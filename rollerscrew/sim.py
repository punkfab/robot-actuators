"""
MuJoCo bench for the roller-screw actuator: one leg of a 57 kg robot, the actuator in its
knee linkage, every number taken from design.py (and the visuals from cad.py's STLs).

What is in the model
  hip       a mass on a vertical slide (half the robot, plus half the payload)
  thigh     rigid with the hip; knee hinge; shank and a foot on a padded, frictionless floor
  actuator  pinned to the thigh at Linkage.base; its rod pinned to the shank's crank
            (a site-to-site connect closes the loop)
     mount ─ series slide (spring IMPACT_K: bearings + housing + load cell) ─ housing
             ─ drive slide (armature = rotor + nut inertia seen through the lead) ─ rod
  screw     forward / back-driving efficiency and seal drag as a load-dependent friction
            budget on the drive slide (friction/ FrictionUpdater, M5 model)
  motor     PID on the rod position from the rotor encoder -> current, clipped by the peak
            current and by the 48 V back-EMF limit, torque capped at the saturation torque;
            winding temperature from copper loss through the design's lumped resistance
  felt      the force in the series spring: what the screw carries and the load cell reads

Experiments (python rollerscrew/sim.py --demo, or make rollerscrew-sim-demo)
  1. statics    hold at four knee angles; spring force against hand statics
  2. landing    drop the leg at 0.5, 1.0, 1.5 m/s with the knee held at 30 deg, for the
                chosen design, the lightest 2 mm design, 3 and 5 mm leads, and no rotor
                inertia; peak screw force against the 12 kN static limit and against
                design.py's one-mass estimate
  3. squat      20 kg payload, 20 -> 100 -> 20 deg in 3 s: force, current, temperature

Modes
  (default)     interactive viewer: sliders for the knee target and the payload
  --demo        the three experiments -> out/sim.png, out/sim.json, out/landing.gif
  --selftest    build, settle, drop once, assert finite and bounded; no window

    make rollerscrew-sim / rollerscrew-sim-demo / rollerscrew-sim-selftest
"""

import os
import sys

_MODE = "demo" if "--demo" in sys.argv else "selftest" if "--selftest" in sys.argv else "interactive"
os.environ.setdefault("MUJOCO_GL", "glfw" if _MODE == "interactive" else "egl")

import json  # noqa: E402
from dataclasses import dataclass, replace  # noqa: E402
from math import atan2, cos, radians, sin, sqrt, tau  # noqa: E402
from pathlib import Path  # noqa: E402

import mujoco  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "friction"))
import design as dz  # noqa: E402
from models import Friction  # noqa: E402
from mujoco_friction import FrictionUpdater  # noqa: E402

OUT = HERE / "out"
CAD = OUT / "cad"
G = 9.81

# leg (assumed; the linkage base sits 180 mm up the thigh)
L_THIGH, L_SHANK, R_FOOT = 0.40, 0.40, 0.03
M_THIGH, M_SHANK = 2.0, 1.2                      # kg, segments without the actuator
HIP_MASS = dz.ROBOT_MASS / 2 - M_THIGH - M_SHANK  # one leg carries half the robot; the actuator adds its own
Y_BONE = 0.06                                     # thigh and shank drawn beside the actuator's plane (visual only)


@dataclass
class Params:
    payload: float = 0.0          # kg on this leg (half the carried load)
    k_series: float = dz.IMPACT_K  # N/m, at the rod
    zeta_series: float = 0.05
    pad_k: float = 3e5            # N/m, foot pad (assumed)
    pad_zeta: float = 0.3
    dt: float = 2e-5
    kp: float = 7.5e5             # N/m on the rod: ~300 N m/rad at the knee with a 20 mm arm
    zeta: float = 0.7
    ki_per_s: float = 4.0         # integral rate, x kp
    ff_gravity: bool = True       # feedforward of the static rod force at the target angle (leg model)
    ff_eta: bool = True           # ... divided by η forward when lifting, times η back when lowering
    ff_load: float = 0.0          # load-cell feedforward gain; 0.9 went unstable in the squat (see README)
    lc_tau: float = 1e-3          # load-cell signal low-pass, s
    no_rotor: bool = False        # A/B: the rotor's inertia removed
    v_bus: float = 48.0


# ---- the actuator, as numbers the model needs ---------------------------------------------
@dataclass
class Act:
    a: dz.Actuator
    label: str

    def __post_init__(self):
        s, m = self.a.screw, self.a.motor
        self.lead = s.lead * 1e-3
        self.k_lead = tau / self.lead                       # rotor rad per m of rod
        self.m_ref, self.j_rot = self.a.reflected_mass()
        mass, parts = self.a.mass()
        self.m_rod = parts["screw"] + 0.02                  # + rod eye
        self.m_housing = mass - self.m_rod
        ph = m.phase()
        self.kt, self.r, self.i_per_j = ph["kt"], ph["r"], ph["i_per_j"]
        self.i_max = self.i_per_j * m.j_peak
        self.t_sat = m.torque(1e9)                          # the saturation cap
        self.eta_f, self.eta_b = s.eta()
        self.drag = s.drag * self.k_lead                    # N at the rod
        self.c_th, self.r_th, self.t_amb = m.heat_capacity(), m.r_th, m.t_amb
        self.od = m.od * 1e-3
        z = self.a.stations()
        self.body_len = (z["front"] - z["cap_rear"]) * 1e-3

    def force_limit(self, i_cmd, v):
        """rod force (N) for a current command (A rms per phase, signed), rod speed v:
        peak-current clip, then the back-EMF limit, then saturation"""
        w = v * self.k_lead
        v_ph = 0.95 * Params.v_bus / sqrt(6)
        e = self.kt * abs(w) / 3
        i = float(np.clip(i_cmd, -self.i_max, self.i_max))
        motoring = i * w > 0
        i_v = max(0.0, (v_ph - e) / self.r) if motoring else (v_ph + e) / self.r
        i = float(np.clip(i, -i_v, i_v))
        t = float(np.clip(self.kt * i, -self.t_sat, self.t_sat))
        return t * self.k_lead, i


def variants():
    """the chosen design and the comparisons, rebuilt by design.py's own sweep"""
    best = dz.chosen()
    req = dz.requirements(best.link)
    rows = dz.sweep(req, best.link)
    ok = lambda r: r["impact"] <= dz.F_PEAK * dz.SHOCK
    by = lambda lead, pred: min((r for r in rows if r["lead"] == lead and pred(r)), key=lambda r: r["mass"])["act"]
    out = [Act(best, "chosen: 2 mm lead, Ø63 motor"),
           Act(by(2.0, lambda r: True), "lightest: 2 mm lead, wider motor"),
           Act(by(3.0, ok), "3 mm lead"),
           Act(by(5.0, ok), "5 mm lead")]
    return out


# ---- geometry ------------------------------------------------------------------------------
def lk_xz(pt):
    """Linkage (x, y) mm, knee origin, thigh +y -> MuJoCo (x, z) m in the knee frame"""
    return pt[0] * 1e-3, pt[1] * 1e-3


def act_dir(link, flex):
    bx, bz = lk_xz(link.base)
    px, pz = lk_xz(link.pin(flex))
    return px - bx, pz - bz


def pose(link, flex_deg, foot_gap=0.0):
    """qpos for knee flexion flex_deg with the foot foot_gap above the floor: hip height,
    actuator pivot angle and rod travel that close the linkage (series spring at rest)"""
    dx0, dz0 = act_dir(link, 0.0)
    dx, dzz = act_dir(link, flex_deg)
    # a positive rotation about +y lowers atan2(z, x)
    pivot = -(atan2(dzz, dx) - atan2(dz0, dx0))
    drive = (link.length(flex_deg) - link.length(0.0)) * 1e-3
    th = radians(flex_deg)
    hip = L_THIGH + L_SHANK * cos(th) + R_FOOT + foot_gap
    return dict(hip=hip, knee=th, pivot=pivot, series=0.0, drive=drive)


MESH_HOUSING = ["housing", "cap_front", "cap_rear", "tail_tube", "stator_core", "coil", "rotor", "magnet", "nut",
                "ring_gear", "carrier", "roller", "bearing_front", "bearing_rear", "encoder_board", "encoder_magnet",
                "load_cell", "clevis"]
MESH_ROD = ["screw", "rod_eye"]
RGBA = dict(housing=".72 .74 .77 1", cap_front=".2 .21 .23 1", cap_rear=".2 .21 .23 1", tail_tube=".2 .21 .23 1",
            stator_core=".25 .27 .3 1", coil=".85 .45 .22 1", rotor=".2 .2 .22 1", magnet=".75 .78 .8 1",
            nut=".58 .52 .42 1", ring_gear=".7 .62 .4 1", carrier=".5 .52 .55 1", roller=".8 .7 .45 1",
            bearing_front=".85 .86 .88 1", bearing_rear=".85 .86 .88 1", encoder_board=".12 .4 .28 1",
            encoder_magnet=".3 .45 .7 1", load_cell=".8 .62 .18 1", clevis=".62 .64 .67 1", screw=".72 .74 .77 1",
            rod_eye=".62 .64 .67 1")


def make_xml(act: Act, p: Params, visual=True):
    link = act.a.link
    bx, bz = lk_xz(link.base)
    cx, cz = lk_xz(link.pin(0.0))
    dx, dzz = act_dir(link, 0.0)
    l0 = link.length(0.0)                                   # mm, pin to pin at flex 0
    pin_min = link.stroke()[1]
    stroke = link.stroke()[0] * 1e-3
    m_ref = 1e-6 if p.no_rotor else act.m_ref
    c_series = 2 * p.zeta_series * sqrt(p.k_series * act.m_rod)
    m_hip = HIP_MASS + p.payload
    tc = 2 * p.dt
    meshes = MESH_HOUSING + MESH_ROD if visual else []
    assets = "\n".join(f'<mesh name="{n}" file="{CAD / (n + ".stl")}" scale=".001 .001 .001"/>' for n in meshes)
    vis = lambda names, z: "\n".join(
        f'<geom type="mesh" mesh="{n}" pos="0 0 {z:.6f}" euler="0 0 90" rgba="{RGBA[n]}" contype="0" conaffinity="0" group="1"/>'
        for n in names if n in meshes)
    r_h = act.od / 2
    zc = act.body_len / 2 + 0.03
    i_rad = act.m_housing * (3 * r_h ** 2 + act.body_len ** 2) / 12
    return f"""<mujoco model="rollerscrew_leg">
  <compiler angle="radian" autolimits="true"/>
  <option timestep="{p.dt}" integrator="implicitfast" gravity="0 0 -{G}"/>
  <visual><global offwidth="1280" offheight="960"/><quality shadowsize="4096"/></visual>
  <asset>
    <texture name="grid" type="2d" builtin="checker" rgb1=".86 .87 .89" rgb2=".78 .8 .82" width="512" height="512"/>
    <material name="grid" texture="grid" texrepeat="8 8"/>
    {assets}
  </asset>
  <worldbody>
    <light pos=".3 -.6 1.6" dir="-.2 .5 -1" diffuse=".9 .9 .9"/>
    <geom name="floor" type="plane" size="2 2 .1" material="grid" contype="0" conaffinity="0"/>
    <body name="hip" pos="0 0 0">
      <joint name="hip" type="slide" axis="0 0 1"/>
      <geom type="box" size=".09 .07 .05" pos="0 0 .05" mass="{m_hip}" rgba=".3 .32 .36 1" contype="0" conaffinity="0"/>
      <geom type="capsule" fromto="0 {Y_BONE} 0 0 {Y_BONE} -{L_THIGH}" size=".03" mass="{M_THIGH}" rgba=".55 .58 .62 1" contype="0" conaffinity="0"/>
      <body name="knee_frame" pos="0 0 -{L_THIGH}">
        <body name="shank">
          <joint name="knee" type="hinge" axis="0 1 0" range="0 {radians(link.flex_max)}"/>
          <geom type="capsule" fromto="0 {Y_BONE} 0 0 {Y_BONE} -{L_SHANK}" size=".025" mass="{M_SHANK}" rgba=".55 .58 .62 1" contype="0" conaffinity="0"/>
          <geom type="cylinder" fromto="0 -.02 0 0 {Y_BONE + 0.03} 0" size=".018" mass="0" rgba=".35 .37 .4 1" contype="0" conaffinity="0"/>
          <geom name="foot" type="sphere" pos="0 0 -{L_SHANK}" size="{R_FOOT}" mass=".01" contype="0" conaffinity="0" rgba=".15 .15 .15 1"/>
          <site name="crank" pos="{cx} 0 {cz}" size=".004" rgba="1 .2 .2 1"/>
          <geom type="capsule" fromto="0 0 0 {cx} 0 {cz}" size=".006" mass="0" rgba=".4 .4 .45 1" contype="0" conaffinity="0"/>
        </body>
        <body name="mount" pos="{bx} 0 {bz}" zaxis="{dx} 0 {dzz}">
          <joint name="pivot" type="hinge" axis="0 1 0"/>
          <inertial pos="0 0 0" mass="0.001" diaginertia="1e-7 1e-7 1e-7"/>
          <body name="housing">
            <joint name="series" type="slide" axis="0 0 1" stiffness="{p.k_series}" damping="{c_series}"/>
            <inertial pos="0 0 {zc}" mass="{act.m_housing}" diaginertia="{i_rad} {i_rad} {act.m_housing * r_h ** 2 / 2}"/>
            {vis(MESH_HOUSING, 0.0)}
            <body name="rod">
              <joint name="drive" type="slide" axis="0 0 1" armature="{m_ref}" range="{-stroke - 0.002} 0.002"/>
              <inertial pos="0 0 {l0 * 1e-3 * 0.75}" mass="{act.m_rod}" diaginertia="1e-5 1e-5 1e-6"/>
              {vis(MESH_ROD, (l0 - pin_min) * 1e-3)}
              <site name="eye" pos="0 0 {l0 * 1e-3}" size=".004" rgba=".2 .4 1 1"/>
            </body>
          </body>
        </body>
      </body>
    </body>
  </worldbody>
  <equality>
    <connect name="loop" site1="eye" site2="crank" solref="{tc} 1"/>
  </equality>
  <actuator>
    <motor name="screw" joint="drive" gear="1" ctrllimited="false"/>
    <motor name="knee_target_deg" joint="knee" gear="0" ctrlrange="0 120"/>
    <motor name="payload_kg" joint="knee" gear="0" ctrlrange="0 40"/>
  </actuator>
</mujoco>"""


class Leg:
    """model + data + controller + friction + winding temperature"""

    def __init__(self, act: Act, p: Params, visual=False):
        self.act, self.p = act, p
        self.m = mujoco.MjModel.from_xml_string(make_xml(act, p, visual))
        self.d = mujoco.MjData(self.m)
        j = lambda n: self.m.joint(n)
        self.q = {n: self.m.jnt_qposadr[j(n).id] for n in ("hip", "knee", "pivot", "series", "drive")}
        self.v = {n: self.m.jnt_dofadr[j(n).id] for n in self.q}
        fr = Friction.from_efficiency(act.eta_f, act.drag, eta_back=act.eta_b)
        self.fr = FrictionUpdater(self.m, "drive", fr)
        self.link = act.a.link
        self.hip_id = self.m.body("hip").id
        self.hip_mass0 = self.m.body_mass[self.hip_id] - p.payload
        self.integ, self.temp, self.i, self.f_lc = 0.0, act.t_amb, 0.0, 0.0
        self.target = 0.0
        self.shank_id = self.m.body("shank").id
        self.foot_gid = self.m.geom("foot").id
        m_leg = HIP_MASS + p.payload + M_THIGH + M_SHANK + act.m_housing + act.m_rod
        self.pad_b = 2 * p.pad_zeta * sqrt(p.pad_k * m_leg)
        self.pad_f = 0.0

    def set_pose(self, flex_deg, foot_gap=0.0):
        for k, v in pose(self.link, flex_deg, foot_gap).items():
            self.d.qpos[self.q[k]] = v
        self.d.qvel[:] = 0
        self.target = self.drive_target(flex_deg)
        self.target_prev = self.target
        self.stance = 1.0 if foot_gap <= 0 else 0.0
        self.integ = 0.0
        self.fr.prime(self.d)
        self.f_lc = self.screw_force()

    def drive_target(self, flex_deg):
        self.target_deg = flex_deg
        return (self.link.length(flex_deg) - self.link.length(0.0)) * 1e-3

    def static_force(self, flex_deg):
        """rod force that holds the leg still at flex_deg, by hand statics on the shank
        (frictionless floor: the foot force is the weight of everything above it)"""
        th = radians(flex_deg)
        above = self.m.body_mass[self.hip_id] + self.act.m_housing + self.act.m_rod     # hip body includes the thigh
        n = (above + M_SHANK + 0.01) * G
        tau = n * L_SHANK * sin(th) - M_SHANK * G * L_SHANK / 2 * sin(th)
        return tau / max(abs(self.link.arm(flex_deg)) * 1e-3, 1e-3)

    def control(self):
        """PID on the rotor encoder (rod travel) plus load-cell feedforward -> current ->
        limited force. A rod pushed in (compression, positive) needs the motor to push out."""
        a, p, d = self.act, self.p, self.d
        x, v = d.qpos[self.q["drive"]], d.qvel[self.v["drive"]]
        m_eff = a.m_ref + a.m_rod if not p.no_rotor else a.m_rod + 1.0
        kd = 2 * p.zeta * sqrt(p.kp * m_eff)
        err = self.target - x
        self.f_lc += (self.screw_force() - self.f_lc) * p.dt / p.lc_tau
        f_cmd = p.kp * err - kd * v + p.kp * p.ki_per_s * self.integ + p.ff_load * self.f_lc
        # stance gate: the gravity feedforward only applies with the foot loaded (a foot switch),
        # ramped over 5 ms so touchdown does not step the command
        self.stance = float(np.clip(self.stance + (1 if self.pad_f > 0 else -1) * p.dt / 5e-3, 0.0, 1.0))
        if p.ff_gravity:
            f_ff = self.static_force(self.target_deg) * self.stance
            if p.ff_eta:
                # target rod speed > 0 extends the rod against the compressive load: driving
                v_t = (self.target - self.target_prev) / p.dt
                w = 0.5 + 0.5 * np.tanh(v_t / 2e-3)              # 0 lowering .. 1 lifting, 2 mm/s blend
                f_ff *= (1 - w) * a.eta_b + w / a.eta_f
            f_cmd += f_ff
        self.target_prev = self.target
        i_cmd = f_cmd / a.k_lead / a.kt
        f, i = a.force_limit(i_cmd, v)
        if abs(i) < a.i_max * 0.999:                        # anti-windup: integrate only when not clipped
            self.integ += err * p.dt
        self.i = i
        d.ctrl[0] = f
        p_cu = 3 * i * i * a.r
        self.temp += (p_cu - (self.temp - a.t_amb) / a.r_th) / a.c_th * p.dt
        return f

    def pad(self):
        """foot pad on a frictionless floor: k x + c xdot while in contact, never pulling.
        Applied here rather than as a MuJoCo contact so k is exactly the stiffness swept."""
        d = self.d
        foot = d.geom_xpos[self.foot_gid]
        pen = R_FOOT - foot[2]
        vel = np.zeros(6)
        mujoco.mj_objectVelocity(self.m, d, mujoco.mjtObj.mjOBJ_GEOM, self.foot_gid, vel, 0)
        f = max(0.0, self.p.pad_k * pen - self.pad_b * vel[5]) if pen > 0 else 0.0
        self.pad_f = f
        d.xfrc_applied[self.shank_id, :] = 0.0
        if f > 0:
            r = foot - d.xipos[self.shank_id]
            d.xfrc_applied[self.shank_id, 2] = f
            d.xfrc_applied[self.shank_id, 3:] = np.cross(r, (0.0, 0.0, f))

    def step(self):
        self.pad()
        self.control()
        self.fr.update(self.d)
        mujoco.mj_step(self.m, self.d)

    # measured quantities
    def screw_force(self):
        """series spring force, compression positive: a rod pushed in shortens the spring
        (q_series < 0), and qfrc_spring = -k q is then positive"""
        return self.d.qfrc_spring[self.v["series"]]

    def flex(self):
        return np.degrees(self.d.qpos[self.q["knee"]])

    def foot_force(self):
        return self.pad_f

    def loop_error(self):
        return float(np.linalg.norm(self.d.site("eye").xpos - self.d.site("crank").xpos))


# ---- experiments ---------------------------------------------------------------------------
def statics(act, p):
    """hold at four angles; screw force against hand statics on the shank"""
    rows = []
    for flex in (20, 45, 70, 90):
        leg = Leg(act, p)
        leg.set_pose(flex)
        for _ in range(int(1.5 / p.dt)):
            leg.step()
        th = radians(leg.flex())
        n = leg.foot_force()
        tau_knee = n * L_SHANK * sin(th) - M_SHANK * G * L_SHANK / 2 * sin(th)
        arm = abs(leg.link.arm(leg.flex())) * 1e-3
        rows.append(dict(flex=leg.flex(), foot=n, sim=leg.screw_force(), hand=tau_knee / arm, i=leg.i,
                         loop_mm=leg.loop_error() * 1e3))
    return rows


def landing(act, p, v_impact, flex=30.0, t_end=0.35, log=False, frames=None):
    """drop from the height that gives v_impact at touchdown; knee held at flex"""
    leg = Leg(act, p, visual=frames is not None)
    h = v_impact ** 2 / (2 * G)
    leg.set_pose(flex, foot_gap=h)
    n = int(t_end / p.dt)
    t_touch = sqrt(2 * h / G) if h > 0 else 0.0
    out = dict(t=[], f=[], flex=[], i=[], v_rod=[], foot=[])
    peak, every = 0.0, max(1, int(1 / 30 / 8 / p.dt))     # gif at 30 fps, 8x slow motion
    renderer = None
    if frames is not None:
        renderer = mujoco.Renderer(leg.m, 720, 960)
        cam = mujoco.MjvCamera()
        knee_z = pose(leg.link, flex)["hip"] - L_THIGH - h
        cam.lookat[:] = (-0.06, 0.0, knee_z - 0.02)
        cam.distance, cam.azimuth, cam.elevation = 0.95, 90, -4
        opt = mujoco.MjvOption()
        opt.geomgroup[:] = 1
    for k in range(n):
        leg.step()
        f = leg.screw_force()
        peak = max(peak, f)
        if log and k % 10 == 0:
            out["t"].append(leg.d.time - t_touch)
            out["f"].append(f)
            out["flex"].append(leg.flex())
            out["i"].append(leg.i)
            out["v_rod"].append(leg.d.qvel[leg.v["drive"]])
            out["foot"].append(leg.foot_force())
        if renderer is not None and k % every == 0 and leg.d.time > t_touch - 0.03:
            renderer.update_scene(leg.d, cam, opt)
            frames.append((renderer.render(), leg.d.time - t_touch, f, leg.foot_force()))
    rod_dv = max(abs(v) for v in out["v_rod"]) if log else None
    return dict(peak=peak, flex_min=None, rod_dv=rod_dv, trace=out if log else None,
                flex_max=max(out["flex"]) if log else None)


def squat(act, p, t_move=3.0, lo=20.0, hi=100.0):
    leg = Leg(act, p)
    leg.set_pose(lo)
    for _ in range(int(0.5 / p.dt)):                       # settle
        leg.step()
    leg.temp = act.t_amb
    t0 = leg.d.time
    out = dict(t=[], target=[], flex=[], f=[], i=[], temp=[], p_el=[])
    for k in range(int((t_move + 0.5) / p.dt)):
        s = min(1.0, (leg.d.time - t0) / t_move)
        u = 0.5 - 0.5 * cos(tau * s)                         # down and up, smooth
        tgt = lo + (hi - lo) * u
        leg.target = leg.drive_target(tgt)
        leg.step()
        if k % 50 == 0:
            out["t"].append(leg.d.time - t0)
            out["target"].append(tgt)
            out["flex"].append(leg.flex())
            out["f"].append(leg.screw_force())
            out["i"].append(leg.i)
            out["temp"].append(leg.temp)
            out["p_el"].append(3 * leg.i ** 2 * act.r + leg.d.ctrl[0] * leg.d.qvel[leg.v["drive"]])
    track = np.abs(np.array(out["flex"]) - np.array(out["target"]))
    return dict(trace=out, f_max=float(np.max(out["f"])), i_max=float(np.max(np.abs(out["i"]))),
                temp_rise=float(np.max(out["temp"]) - act.t_amb), track_max=float(track.max()),
                track_rms=float(np.sqrt(np.mean(track ** 2))))


def demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    OUT.mkdir(exist_ok=True)
    acts = variants()
    p = Params()
    chosen = acts[0]
    print(f"chosen: housing {chosen.m_housing:.3f} kg, rod {chosen.m_rod:.3f} kg, rotor seen at the rod {chosen.m_ref:.0f} kg, "
          f"peak {chosen.i_max:.1f} A, saturation {chosen.t_sat:.2f} N m = {chosen.t_sat * chosen.k_lead:.0f} N")

    print("\n1. statics: screw force held at each angle, against hand statics on the shank")
    st = statics(chosen, p)
    for r in st:
        print(f"   knee {r['flex']:5.1f} deg  foot {r['foot']:6.0f} N  screw sim {r['sim']:7.1f} N  hand {r['hand']:7.1f} N "
              f"({r['sim'] / r['hand'] - 1:+.2%})  current {r['i']:5.2f} A  loop gap {r['loop_mm']:.4f} mm")

    print("\n2. landing, knee held at 30 deg (peak screw force, kN; static limit "
          f"{dz.F_PEAK * dz.SHOCK / 1e3:.0f} kN)")
    speeds = (0.5, 1.0, 1.5)
    cases = [(a, replace(p), a.label) for a in acts] + [(chosen, replace(p, no_rotor=True), "chosen, no rotor inertia"),
                                                        (chosen, replace(p, ff_gravity=False), "chosen, PID only (no feedforward)")]
    land = []
    for a, pp, label in cases:
        row = dict(label=label, m_ref=0.0 if pp.no_rotor else a.m_ref, peaks=[], estimate=dz.landing(a.a) if not pp.no_rotor else None)
        for v in speeds:
            r = landing(a, pp, v, log=True)
            row["peaks"].append(r["peak"])
            row.setdefault("rod_dv", []).append(r["rod_dv"])
            if v == 1.0:
                row["trace"] = r["trace"]
        land.append(row)
        print(f"   {label:34} rotor {row['m_ref']:4.0f} kg: " + "  ".join(
            f"{v} m/s {pk / 1e3:5.1f}" for v, pk in zip(speeds, row["peaks"]))
              + f"   rod speed change at 1 m/s {row['rod_dv'][1]:.2f} m/s"
              + (f"   (design.py one-mass {row['estimate'] / 1e3:.1f} kN at {dz.IMPACT_V} m/s)" if row["estimate"] else ""))
    soft = landing(chosen, replace(p, pad_k=1e5), 1.0)["peak"]
    hard = landing(chosen, replace(p, pad_k=1e6), 1.0)["peak"]
    print(f"   foot pad A/B, chosen at 1 m/s: soft 1e5 N/m {soft / 1e3:.1f} kN, as run 3e5 {land[0]['peaks'][1] / 1e3:.1f} kN, "
          f"hard 1e6 {hard / 1e3:.1f} kN")

    print("\n3. squat, 20 kg payload (10 kg on this leg), 20 -> 100 -> 20 deg in 3 s")
    ctrls = [("PID only", dict(ff_gravity=False)), ("+ gravity feedforward", dict(ff_eta=False)),
             ("+ gravity / screw efficiency", dict())]
    sqs = {}
    for name, kw in ctrls:
        r = squat(chosen, replace(p, payload=dz.PAYLOAD / 2, **kw))
        sqs[name] = r
        print(f"   {name:30} peak screw {r['f_max']:5.0f} N, peak current {r['i_max']:4.1f} A of {chosen.i_max:.1f}, "
              f"winding rise {r['temp_rise']:.2f} K, tracking max {r['track_max']:5.2f} deg rms {r['track_rms']:.2f} deg")
    sq = sqs[ctrls[-1][0]]
    probe = Leg(chosen, replace(p, payload=dz.PAYLOAD / 2))
    print(f"   hand statics in this leg at 100 deg: {probe.static_force(100):.0f} N (ground force on a "
          f"{L_SHANK * sin(radians(100)):.2f} m lever, thigh upright); design.py's duty assumed a 0.12 m lever: "
          f"{dz.requirements(chosen.a.link)['f_squat']:.0f} N with its x1.5 rising factor")

    # figure
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, ax = plt.subplots(2, 2, figsize=(11, 7.6))
    cols = ["#1f4e79", "#b03a2e", "#2a7f3f", "#6a51a3", "0.55", "#c87a2f"]
    a0 = ax[0, 0]
    for row, c in zip(land, cols):
        tr = row["trace"]
        a0.plot(np.array(tr["t"]) * 1e3, np.array(tr["f"]) / 1e3, color=c, lw=1.3, label=row["label"])
    a0.axhline(dz.F_PEAK * dz.SHOCK / 1e3, color="k", lw=0.7, ls="--")
    a0.text(2, dz.F_PEAK * dz.SHOCK / 1e3 + 0.2, "screw static limit", fontsize=8)
    a0.set(xlim=(-5, 80), xlabel="time from touchdown (ms)", ylabel="screw force (kN)", title="landing at 1 m/s, knee held at 30°")
    a0.legend(frameon=False, fontsize=7.5, loc="upper right", bbox_to_anchor=(1.0, 0.93))
    a1 = ax[0, 1]
    for row, c in zip(land, cols):
        a1.plot(speeds, np.array(row["peaks"]) / 1e3, "o-", color=c, ms=4, lw=1.3, label=row["label"])
    a1.axhline(dz.F_PEAK * dz.SHOCK / 1e3, color="k", lw=0.7, ls="--")
    a1.set(xlabel="touchdown speed (m/s)", ylabel="peak screw force (kN)", title="peak screw force against touchdown speed")
    a2 = ax[1, 0]
    tr = sq["trace"]
    a2.plot(tr["t"], tr["target"], color="0.6", lw=1.2, ls="--", label="target")
    t0 = sqs[ctrls[0][0]]["trace"]
    a2.plot(t0["t"], t0["flex"], color="#1f4e79", lw=1.0, alpha=0.45, label="knee, PID only")
    a2.plot(tr["t"], tr["flex"], color="#1f4e79", lw=1.4, label="knee, with feedforward")
    a2.set(xlabel="time (s)", ylabel="knee flexion (deg)", title="squat with 20 kg (10 kg on this leg)")
    q = a2.twinx()
    q.spines["right"].set_visible(True)
    q.plot(tr["t"], np.array(tr["f"]) / 1e3, color="#b03a2e", lw=1.2, label="screw force")
    q.set_ylabel("screw force (kN)", color="#b03a2e")
    a2.legend(frameon=False, fontsize=8, loc="upper left")
    a3 = ax[1, 1]
    a3.plot(tr["t"], tr["i"], color="#c87a2f", lw=1.3, label="phase current (A rms)")
    a3.axhline(chosen.i_max, color="#c87a2f", lw=0.7, ls="--")
    a3.set(xlabel="time (s)", ylabel="current (A)", title="squat: current and winding temperature")
    q = a3.twinx()
    q.spines["right"].set_visible(True)
    q.plot(tr["t"], np.array(tr["temp"]) - chosen.t_amb, color="#6a51a3", lw=1.3)
    q.set_ylabel("winding rise (K)", color="#6a51a3")
    fig.tight_layout()
    fig.savefig(OUT / "sim.png", dpi=150)

    frames = []
    landing(chosen, p, 1.0, t_end=0.25, frames=frames)
    from PIL import Image, ImageDraw
    fr = []
    for img, t, f, foot in frames:
        im = Image.fromarray(img)
        dr = ImageDraw.Draw(im)
        bar = int(np.clip(f / (dz.F_PEAK * dz.SHOCK), 0, 1.2) * 600)
        dr.rectangle((40, 660, 40 + 600, 690), outline=(200, 200, 200))
        dr.rectangle((40, 660, 40 + bar, 690), fill=(200, 60, 50) if f > dz.F_PEAK * dz.SHOCK else (230, 150, 60))
        dr.text((40, 620), f"t = {t * 1e3:6.1f} ms from touchdown   screw {f / 1e3:5.1f} kN   foot {foot / 1e3:4.2f} kN"
                f"   (bar full = {dz.F_PEAK * dz.SHOCK / 1e3:.0f} kN static limit)", fill=(240, 240, 240))
        dr.text((40, 30), f"{chosen.label}: 1 m/s landing, knee held at 30 deg, 8x slow motion", fill=(240, 240, 240))
        fr.append(im.convert("P", palette=Image.Palette.ADAPTIVE, colors=96))
    fr[0].save(OUT / "landing.gif", save_all=True, append_images=fr[1:], duration=33, loop=0)
    res = dict(statics=st, landing=[{k: v for k, v in r.items() if k != "trace"} for r in land], speeds=speeds,
               pad_ab=dict(soft=soft, nominal=land[0]["peaks"][1], hard=hard),
               squat={n: {k: v for k, v in r.items() if k != "trace"} for n, r in sqs.items()}, params=vars(p))
    (OUT / "sim.json").write_text(json.dumps(res, indent=1, default=float))
    print(f"\nwrote {OUT / 'sim.png'}, {OUT / 'landing.gif'} ({len(frames)} frames), {OUT / 'sim.json'}")


def selftest():
    act = Act(dz.chosen(), "chosen")
    p = Params()
    leg = Leg(act, p)
    leg.set_pose(30)
    print(f"model: {leg.m.nbody} bodies, {leg.m.njnt} joints, masses (kg): " + ", ".join(
        f"{leg.m.body(i).name} {leg.m.body_mass[i]:.3f}" for i in range(1, leg.m.nbody)))
    print(f"armature on the drive slide {leg.m.dof_armature[leg.v['drive']]:.0f} kg; loop gap at start {leg.loop_error()*1e3:.4f} mm")
    for _ in range(int(0.5 / p.dt)):
        leg.step()
    assert np.all(np.isfinite(leg.d.qpos)), "non-finite state"
    assert leg.loop_error() < 1e-4, f"linkage loop opened: {leg.loop_error()}"
    assert 25 < leg.flex() < 35, f"knee drifted to {leg.flex():.1f}"
    r = landing(act, p, 1.0, t_end=0.2)
    assert np.isfinite(r["peak"]) and 0 < r["peak"] < 1e5, r["peak"]
    # the interactive model: CAD meshes and the coarser real-time step
    pv = Params(dt=1e-4)
    lv = Leg(act, pv, visual=True)
    lv.set_pose(45)
    for _ in range(int(1.0 / pv.dt)):
        lv.step()
    assert np.all(np.isfinite(lv.d.qpos)) and abs(lv.flex() - 45) < 3 and lv.loop_error() < 1e-3, \
        (lv.flex(), lv.loop_error())
    print(f"interactive model ({lv.m.nmesh} CAD meshes, dt {pv.dt}): held 45 deg at {lv.flex():.2f}")
    print(f"selftest ok: held 30 deg at {leg.flex():.2f}, screw {leg.screw_force():.0f} N; 1 m/s landing peak {r['peak']/1e3:.1f} kN")


def interactive():
    import time

    import mujoco.viewer
    act = Act(dz.chosen(), "chosen")
    p = Params(dt=1e-4)                                    # coarser for real time; the batch runs use 2e-5
    leg = Leg(act, p, visual=True)
    leg.set_pose(30)
    m, d = leg.m, leg.d
    i_t, i_p = m.actuator("knee_target_deg").id, m.actuator("payload_kg").id
    d.ctrl[i_t], d.ctrl[i_p] = 30.0, 0.0
    print("sliders: knee_target_deg, payload_kg. Ctrl+right-drag on the hip to push it. "
          "Screw force and winding temperature print once a second.")
    with mujoco.viewer.launch_passive(m, d) as v:
        last = time.time()
        while v.is_running():
            t0 = time.time()
            m.body_mass[leg.hip_id] = leg.hip_mass0 + d.ctrl[i_p]
            leg.target = leg.drive_target(float(d.ctrl[i_t]))
            tgt, pl = d.ctrl[i_t], d.ctrl[i_p]
            for _ in range(int(1 / 60 / p.dt)):
                leg.step()
                d.ctrl[i_t], d.ctrl[i_p] = tgt, pl
            v.sync()
            if time.time() - last > 1:
                last = time.time()
                print(f"knee {leg.flex():5.1f} deg  screw {leg.screw_force():7.0f} N  current {leg.i:5.1f} A  "
                      f"winding {leg.temp:5.1f} C")
            time.sleep(max(0.0, 1 / 60 - (time.time() - t0)))


if __name__ == "__main__":
    {"interactive": interactive, "demo": demo, "selftest": selftest}[_MODE]()
