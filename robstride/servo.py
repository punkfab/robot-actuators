"""
RobStride actuators in MuJoCo: the driver's MIT law, not an ideal position servo.

A RobStride is a quasi-direct-drive (QDD) joint: an outrunner BLDC, one planetary
stage (~7–10:1) and a FOC driver. The host sends MIT frames over CAN,

    (p*, v*, kp, kd, τff)  ->  τ = kp·(p* − q) + kd·(v* − q̇) + τff

and the driver closes that law at its own rate against its own encoder. Next to a
hobby servo (../sts3215) the law is simple and fully exposed: kp IS the joint
stiffness in N·m/rad, there's no hidden firmware gain. What a MuJoCo
<position kp=…> leaves out is everything around the law:

  command path   host rate (ZOH) + CAN latency; every field quantized to the
                 packet's bits (p 16, v/kp/kd/τ 12) over the model's ranges
  torque limits  ±τ_max of the MIT range, and the published 48 V T-N curve: at
                 speed the bus voltage runs out and τ falls to 0 at no-load speed
                 (speed axis scaled by V_bus)
  gearbox        one planetary stage: load-dependent friction from ../friction,
                 in the form Katz measured on the mini-cheetah actuator, until a
                 bench fits it (bench.py)
  rotor          reflected rotor inertia N²·J_r as armature
  thermal        two-node winding model fitted to the overload tables (thermal.py):
                 θ = 1 is where the driver's protection trips (torque goes to 0
                 here, as the real one faults); holding still heats harder than
                 rotating

Usage (same pattern as sts3215/arm.py):

    motorize(spec, {"hip": "RS03", "knee": "RS03"})     # joints -> torque motors
    m = spec.compile(); d = mujoco.MjData(m)
    rs = RobStrideServos(m, d, {"hip": "RS03", "knee": "RS03"}, kp=80, kd=4)
    rs.set("hip", p=0.4)                  # optional v=, tau_ff=, kp=, kd=
    rs.step()                             # law + friction + mj_step

Numbers for each model live in specs.py, with their sources.
"""

import sys
from dataclasses import dataclass, field
from pathlib import Path

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "friction"))
from models import Friction  # noqa: E402
from mujoco_friction import FrictionUpdater  # noqa: E402
from specs import RPM, SPECS, Spec  # noqa: E402
from thermal import Thermal  # noqa: E402


def quantize(x, lo, hi, bits):
    """Round-trip through the MIT packet's unsigned int field (clips to the range)."""
    n = (1 << bits) - 1
    u = np.round((np.clip(x, lo, hi) - lo) / (hi - lo) * n)
    return lo + u * (hi - lo) / n


def envelope(spec: Spec, dq, tau, vbus=None):
    """Largest |τ| the driver can make at output speed dq, for a command of sign(τ).
    Motoring (τ·dq > 0) follows the booklet's T-N table, with the speed axis scaled by
    V_bus/V_nom (back-EMF ∝ speed); braking is limited by current only (τ_peak)."""
    vbus = spec.v_nom if vbus is None else vbus
    t_tab, n_tab = spec.tn_arrays()
    w = np.abs(dq) / RPM * spec.v_nom / vbus         # equivalent speed at V_nom
    t_v = np.interp(w, n_tab, t_tab, right=0.0)
    motoring = np.sign(tau) * np.sign(dq) > 0
    return np.where(motoring, t_v, spec.tau_peak)


def gearbox_prior(spec: Spec, friction_scale=1.0) -> Friction:
    """Planetary prior at the output, in Katz's measured form (thesis Eq. 2.27):
    budget = drag + k_load·|τ_motor|, a BAM M5 with Ke = 0, so η_drive = 1 − k_load
    and η_back = 1/(1 + k_load), ≈ 0.96 both ways. friction_scale multiplies the whole
    budget (drag and load terms), as a DR knob."""
    f = Friction(kc=spec.drag, kv=spec.kv, km=spec.k_load, ke=0.0)
    f = f.scaled(kc=f.kc * friction_scale, km=f.km * friction_scale,
                 ke=f.ke * friction_scale, kv=f.kv * friction_scale)
    f.name = f"{spec.name} prior"
    return f


def motorize(spec, joints: dict, name_fmt="{j}"):
    """Replace/insert a plain torque motor (gear 1, unlimited ctrl) on each joint.
    Existing actuators named name_fmt.format(j=joint) are converted; others added.
    Joint limits stay whatever the model says (the MIT law has no position limit)."""
    existing = {a.name: a for a in spec.actuators}
    for j in joints:
        a = existing.get(name_fmt.format(j=j))
        if a is None:
            a = spec.add_actuator(name=name_fmt.format(j=j), target=j,
                                  trntype=mujoco.mjtTrn.mjTRN_JOINT)
        a.set_to_motor()
        a.gear[0] = 1.0
        a.ctrllimited = mujoco.mjtLimited.mjLIMITED_FALSE
        a.forcelimited = mujoco.mjtLimited.mjLIMITED_FALSE


@dataclass
class Joint:
    spec: Spec
    dof: int
    qadr: int
    aid: int
    fr: FrictionUpdater
    kp: float
    kd: float
    p: float = 0.0
    v: float = 0.0
    tau_ff: float = 0.0
    applied: tuple = (0.0, 0.0, 0.0, 0.0, 0.0)   # the packet the driver is running
    queue: list = field(default_factory=list)      # (t_arrive, packet)
    thermal: Thermal = None                        # θ: 1 = protection limit
    tau: float = 0.0                               # last output torque (N·m)
    fault: float = None                            # sim time the over-temp trip fired


class RobStrideServos:
    """MIT-mode RobStride drivers on joints of a compiled model (after motorize()).

    joints        {joint_name: "RS03", ...}  (actuator name = joint name by default)
    kp, kd        default MIT gains, N·m/rad and N·m·s/rad, output side
    host_hz       rate the host sends frames (ZOH between them)
    latency       s, host -> driver (CAN + driver pickup)
    vbus          V; the back-EMF envelope scales with it
    friction_scale / armature_scale      real2sim multipliers and DR knobs
    friction      Friction (all joints) or {joint: Friction}: a fitted model
                  (bench.py) instead of the spec prior; armature likewise (kg·m²)
    ideal=True    no quantization, latency, envelope, friction or armature: the
                  "MuJoCo <position>" control case, for comparisons
    """

    def __init__(self, model, data, joints: dict, kp=40.0, kd=2.0, host_hz=500.0,
                 latency=0.5e-3, vbus=None, friction_scale=1.0, armature_scale=1.0,
                 ideal=False, name_fmt="{j}", friction=None, armature=None):
        self.m, self.d = model, data
        self.host_dt = 1.0 / host_hz
        self.latency = latency
        self.vbus = vbus
        self.ideal = ideal
        self.next_send = 0.0
        self.j = {}
        for name, kind in joints.items():
            s = SPECS[kind] if isinstance(kind, str) else kind
            jid = model.joint(name).id
            dof = model.jnt_dofadr[jid]
            if ideal:
                fr = FrictionUpdater(model, name, Friction(), solimp=None)
            else:
                f = friction.get(name) if isinstance(friction, dict) else friction
                fr = FrictionUpdater(model, name, f or gearbox_prior(s, friction_scale))
                a = armature.get(name) if isinstance(armature, dict) else armature
                model.dof_armature[dof] = (a if a is not None else s.armature) * armature_scale
            self.j[name] = Joint(s, dof, model.jnt_qposadr[jid],
                                 model.actuator(name_fmt.format(j=name)).id, fr, kp, kd,
                                 thermal=Thermal(s))
            self.set(name, p=float(data.qpos[model.jnt_qposadr[jid]]))

    # ---- host side ----------------------------------------------------------------
    def set(self, name, p=None, v=None, tau_ff=None, kp=None, kd=None):
        j = self.j[name]
        for k, x in dict(p=p, v=v, tau_ff=tau_ff, kp=kp, kd=kd).items():
            if x is not None:
                setattr(j, k, float(x))

    def _packet(self, j: Joint):
        s = j.spec
        if self.ideal:
            return (j.p, j.v, j.kp, j.kd, j.tau_ff)
        return (quantize(j.p, -s.p_max, s.p_max, 16),
                quantize(j.v, -s.v_max, s.v_max, 12),
                quantize(j.kp, 0.0, s.kp_max, 12),
                quantize(j.kd, 0.0, s.kd_max, 12),
                quantize(j.tau_ff, -s.t_max, s.t_max, 12))

    def reset(self):
        self.next_send = self.d.time
        for j in self.j.values():
            j.queue.clear()
            j.thermal.reset()
            j.fault = None
            j.applied = self._packet(j)

    # ---- driver side --------------------------------------------------------------
    def update(self):
        t = self.d.time
        if t + 1e-12 >= self.next_send:
            self.next_send += self.host_dt
            for j in self.j.values():
                pkt = self._packet(j)
                if self.ideal or self.latency <= 0:
                    j.applied = pkt
                else:
                    j.queue.append((t + self.latency, pkt))
        for j in self.j.values():
            while j.queue and j.queue[0][0] <= t + 1e-12:
                j.applied = j.queue.pop(0)[1]
            p, v, kp, kd, tff = j.applied
            q, dq = self.d.qpos[j.qadr], self.d.qvel[j.dof]
            tau = kp * (p - q) + kd * (v - dq) + tff
            if not self.ideal:
                s = j.spec
                tau = float(np.clip(tau, -s.t_max, s.t_max))
                lim = float(envelope(s, dq, tau, self.vbus))
                tau = float(np.clip(tau, -lim, lim))
                if j.fault is None and j.thermal.step(tau, dq, self.m.opt.timestep) >= 1.0:
                    j.fault = t                    # driver trips: torque off until reset
                if j.fault is not None:
                    tau = 0.0
            j.tau = tau
            self.d.ctrl[j.aid] = tau
            j.fr.update(self.d)

    def step(self, n=1):
        for _ in range(n):
            self.update()
            mujoco.mj_step(self.m, self.d)

    def state(self, name):
        j = self.j[name]
        return dict(q=float(self.d.qpos[j.qadr]), dq=float(self.d.qvel[j.dof]),
                    tau=j.tau, theta=j.thermal.theta)
