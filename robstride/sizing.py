"""
Pick RobStride models for a design's joints from what the motion actually demands.

Datasheet selection ("peak torque > gravity load") misses three things that decide
QDD sizing in practice:

  1. Heat, not peak. The rated torque is the thermal limit; a cycle that sits at
     60% of peak is fine for 2 s and trips the protection in under a minute. And
     holding still is worse than moving (thermal.py): RS03 is rated 21 N·m turning,
     13 N·m holding. The cycle is checked on the fitted thermal model, repeated
     forever. θ scales with torque², so the margin applies squared.
  2. The speed envelope. Near no-load speed the bus voltage, not the current, limits
     torque, so a fast move can fail at half of peak torque.
  3. The actuators are part of the load. Each one weighs 0.3–2.4 kg and sits on the
     link before it, so the choice for the wrist changes the shoulder's demand.
     sizing iterates: choose -> add their masses -> re-run inverse dynamics.

Demand per joint comes from MuJoCo inverse dynamics on the design's own MJCF along a
trajectory, plus the candidate's reflected armature (N²·J_r·q̈) and gearbox friction
prior. Nothing here is a controller: it's the torque the joint must make IF it
follows the trajectory. servo.py checks what the MIT loop actually does.

    from sizing import size_joints
    picks = size_joints(spec_fn, joints, traj, dt)          # see demo.py
"""

import sys
from dataclasses import dataclass
from pathlib import Path

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from servo import envelope, gearbox_prior  # noqa: E402
from specs import SPECS, Spec  # noqa: E402
from thermal import cycle_theta  # noqa: E402


@dataclass
class Check:
    kind: str
    peak: float        # max |τ| / τ_peak
    theta: float       # thermal load, cycle repeated forever (1 = protection trips)
    env: float         # max over the cycle of |τ| / T-N envelope at that speed
    speed: float       # max |ω| / no-load speed
    tau_peak: float    # N·m demanded
    tau_rms: float
    ok: bool

    def row(self):
        flag = "ok " if self.ok else "NO "
        return (f"{flag}{self.kind:5s} peak {self.tau_peak:6.1f} N·m ({self.peak:4.0%})  "
                f"rms {self.tau_rms:5.1f}  thermal {self.theta:4.0%}  "
                f"T-N {self.env:4.0%}  speed {self.speed:4.0%}")


def inverse_dynamics(model, qpos, dt):
    """qfrc_inverse along a joint-space trajectory qpos[T, nq] (hinges only)."""
    qvel = np.gradient(qpos, dt, axis=0)
    qacc = np.gradient(qvel, dt, axis=0)
    d = mujoco.MjData(model)
    out = np.zeros_like(qpos)
    for k in range(len(qpos)):
        d.qpos[:], d.qvel[:], d.qacc[:] = qpos[k], qvel[k], qacc[k]
        mujoco.mj_inverse(model, d)
        out[k] = d.qfrc_inverse
    return out, qvel, qacc


def check(spec: Spec, tau_load, dq, ddq, dt, margin=1.0, vbus=None):
    """tau_load: load torque at the joint (no actuator losses). Adds armature and the
    gearbox prior (motor must cover friction too) and checks every limit."""
    f = gearbox_prior(spec)
    tau = tau_load + spec.armature * ddq
    # Gearbox: driving (τ·ω > 0) the motor pays τ/η_f, backdriven (load helping) it
    # only has to hold τ·η_b; drag and viscous friction always oppose motion. At rest
    # friction helps hold, so the motor needs at most the load itself.
    drive = tau * dq > 0
    tau_m = np.where(drive, tau / f.eta_drive, tau * f.eta_back) \
        + np.sign(dq) * (spec.drag + f.kv * np.abs(dq))
    peak = np.max(np.abs(tau_m))
    rms = float(np.sqrt(np.mean(tau_m ** 2)))
    env = np.max(np.abs(tau_m) / np.maximum(envelope(spec, dq, tau_m, vbus), 1e-9))
    speed = np.max(np.abs(dq)) / (spec.w_noload * (vbus or spec.v_nom) / spec.v_nom)
    theta, theta_peak, _ = cycle_theta(spec, tau_m, dq, dt)
    c = Check(spec.name, peak / spec.tau_peak, max(theta, theta_peak), float(env),
              float(speed), float(peak), rms, False)
    c.ok = c.peak <= 1 / margin and c.theta <= 1 / margin**2 and c.env <= 1 / margin
    return c


def add_actuator_masses(spec, picks: dict, stator_on_parent=True):
    """Put each chosen actuator's mass at its joint, on the parent link (stator side),
    as a cylinder of the model's diameter and length along the joint axis."""
    for jname, kind in picks.items():
        s = SPECS[kind]
        j = next(x for x in spec.joints if x.name == jname)
        body = j.parent.parent if stator_on_parent and j.parent.parent is not None else j.parent
        # joint anchor expressed in the parent: child pos + joint pos (no child rotation assumed)
        pos = np.array(j.parent.pos) + np.array(j.pos)
        axis = np.array(j.axis, float) / np.linalg.norm(j.axis)
        z = np.array([0, 0, 1.0])
        v = np.cross(z, axis)
        quat = [1, 0, 0, 0] if np.linalg.norm(v) < 1e-9 else np.r_[
            np.cos(np.arccos(np.clip(z @ axis, -1, 1)) / 2),
            np.sin(np.arccos(np.clip(z @ axis, -1, 1)) / 2) * v / np.linalg.norm(v)]
        body.add_geom(name=f"rs_{jname}", type=mujoco.mjtGeom.mjGEOM_CYLINDER,
                      size=[s.od / 2, s.length / 2, 0], pos=pos, quat=quat,
                      mass=s.mass, contype=0, conaffinity=0, rgba=[0.15, 0.15, 0.17, 1])


def size_joints(make_spec, joints, qpos, dt, candidates=None, margin=1.2,
                vbus=None, max_iter=6, verbose=True):
    """make_spec() -> fresh MjSpec of the design (without actuator masses).
    joints: joint names to size, in qpos order of the trajectory columns (hinges).
    Returns (picks, checks, every): lightest passing model per joint, its Check, and
    the Checks of every candidate (why the lighter ones failed)."""
    candidates = candidates or sorted(SPECS, key=lambda k: SPECS[k].mass)
    picks = {j: candidates[0] for j in joints}
    for it in range(max_iter):
        spec = make_spec()
        add_actuator_masses(spec, picks)
        m = spec.compile()
        tau, dq, ddq = inverse_dynamics(m, qpos, dt)
        new, checks, every = {}, {}, {}
        for i, jn in enumerate(joints):
            dof = m.jnt_dofadr[m.joint(jn).id]
            every[jn] = [check(SPECS[k], tau[:, dof], dq[:, dof], ddq[:, dof], dt, margin, vbus)
                         for k in candidates]
            c = next((c for c in every[jn] if c.ok), every[jn][-1])
            new[jn], checks[jn] = c.kind, c
        if verbose:
            print(f"  iter {it}: " + "  ".join(f"{j}={k}" for j, k in new.items()))
        if new == picks:
            break
        picks = new
    return picks, checks, every


def check_picks(make_spec, joints, picks, qpos, dt, margin=1.0, vbus=None):
    """Check a given choice (e.g. a datasheet guess) against the same demand."""
    spec = make_spec()
    add_actuator_masses(spec, picks)
    m = spec.compile()
    tau, dq, ddq = inverse_dynamics(m, qpos, dt)
    out = {}
    for jn in joints:
        dof = m.jnt_dofadr[m.joint(jn).id]
        out[jn] = check(SPECS[picks[jn]], tau[:, dof], dq[:, dof], ddq[:, dof], dt, margin, vbus)
    return out
