"""
MuJoCo (4.0) model of wobble2, from the CAD meshes and design.py.

Bodies: the input carrier (hinge, driven by a current- and voltage-limited motor), the wobble ring
(a ball joint at the nutation centre for the spherical seat, on a light axial slide for the seat's
axial clearance), and the output carrier (hinge on its main bearing) with a load arm.

The roller contacts are a custom force law, as in wobble/: MuJoCo would collide a roller with the
wavy face's convex hull. Each step, for the 15 rollers: the approach of roller into face from the
ring's and output's deviation from their conjugate poses at the current input angle (design.py's
contact Jacobian), a Palmgren line-contact force with damping along the face normal, and rolling
friction along the path. The pins hold the ring's rotation about the bisector of the housing and
ring axes (a constant-velocity constraint) and slide with friction. Forces go in between mj_step1
and mj_step2 so they see the state they act on.

    python wobble2/sim.py              interactive viewer (slow motion; slider: output target)
    python wobble2/sim.py --demo       ratio, hold, efficiency, shaking -> out/sim_demo.png, sim.json, sim.gif
    python wobble2/sim.py --selftest   build, hold a load, compare with design.solve
"""

import os
import sys

_MODE = "demo" if "--demo" in sys.argv else "selftest" if "--selftest" in sys.argv else "interactive"
os.environ.setdefault("MUJOCO_GL", "glfw" if _MODE == "interactive" else "egl")

import json  # noqa: E402
from dataclasses import dataclass, replace  # noqa: E402
from pathlib import Path  # noqa: E402

import mujoco  # noqa: E402
import numpy as np  # noqa: E402
import trimesh  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
SIM = OUT / "sim"
sys.path.insert(0, str(HERE))
import cad as C  # noqa: E402
import design as D  # noqa: E402

# world <- CAD: CAD x -> world z (up), CAD y -> world x, CAD z -> world y: the axis lies horizontal
G = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0.0]])
CUTAWAY = False


@dataclass
class Params:
    dt: float = 1e-5
    c_contact: float = 300.0         # N s/m per roller contact
    k_pin: float = 2e5               # N m/rad pins about the bisector
    c_pin: float = 5.0
    k_ax: float = 5e4                # N/m the seat's axial float (design.solve k_ax)
    arm_mass: float = 4.0
    arm_len: float = 0.30
    Kt: float = 0.15                 # N m/A  (stand-in motor: 4.5 N m peak)
    R: float = 0.15
    V: float = 48.0
    I_max: float = 30.0
    J_rotor: float = 3e-4
    kp: float = 1500.0               # N m/rad on the output angle
    kd: float = 40.0
    k_lock: float = 600.0            # N m/rad input held ("lock" mode, the statics check)
    c_lock: float = 2.0
    out_damping: float = 0.0         # N m s/rad on the output hinge (implicit; settles the statics check)
    load_torque: float = 0.0         # N m, a dyno: constant torque against the output's motion


def quat(Rm):
    q = np.zeros(4)
    mujoco.mju_mat2Quat(q, np.asarray(Rm, float).flatten())
    return q


def rotvec(Rm):
    c = np.clip((np.trace(Rm) - 1) / 2, -1, 1)
    a = np.arccos(c)
    v = np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]])
    return v / 2 if a < 1e-9 else v * a / (2 * np.sin(a))


def ring_rate(w, phi, h=1e-6):
    A = D.ring_rot(w, phi + h) @ D.ring_rot(w, phi - h).T
    Wm = (A - A.T) / 2 / (2 * h)
    return np.array([Wm[2, 1], Wm[0, 2], Wm[1, 0]])


def mesh_files():
    SIM.mkdir(parents=True, exist_ok=True)
    cad = OUT / "cad"
    if not (cad / "ring.stl").exists():
        C.main()
    f = {p.stem: str(p) for p in cad.glob("*.stl")}
    if CUTAWAY:
        box = trimesh.creation.box(extents=(400, 400, 400))
        box.apply_translation((0, -200, 0))
        for n in ("housing_in", "housing_out"):
            C.boolean(trimesh.load(f[n]), box, "diff").export(SIM / f"{n}_cut.stl")
            f[n] = str(SIM / f"{n}_cut.stl")
    return f


def make_xml(w, p, files):
    asset = "\n".join(f'<mesh name="{n}" file="{fn}" scale="0.001 0.001 0.001" inertia="exact"/>' for n, fn in files.items())
    gq = " ".join(map(str, quat(G)))

    def rollers(name, row):
        return "".join(f'<geom class="vis" type="mesh" mesh="{name}" quat="{" ".join(map(str, quat(C.rot_to(u))))}" '
                       f'rgba=".85 .86 .9 1" density="7850"/>' for u in D.axes(row))
    oq = quat(D.Rz(D.out_angle(w, 0.0)))
    return f"""<mujoco model="wobble2">
  <compiler angle="radian"/>
  <option timestep="{p.dt}" integrator="Euler" gravity="0 0 -9.81"/>
  <visual><global offwidth="1280" offheight="720"/><headlight ambient=".4 .4 .4"/></visual>
  <default><default class="vis"><geom contype="0" conaffinity="0"/></default></default>
  <asset>{asset}
    <texture name="sky" type="skybox" builtin="gradient" rgb1=".97 .97 .98" rgb2=".78 .79 .82" width="64" height="64"/>
  </asset>
  <worldbody>
    <light pos="0 -0.4 0.6" dir="0 .5 -1" diffuse=".8 .8 .8"/>
    <body name="frame" quat="{gq}">
      <geom class="vis" type="mesh" mesh="housing_in" rgba=".8 .81 .83 .35"/>
      <geom class="vis" type="mesh" mesh="housing_out" rgba=".8 .81 .83 .35"/>
      <geom class="vis" type="mesh" mesh="motor" rgba=".12 .12 .13 1"/>
      <body name="input">
        <joint name="input" type="hinge" axis="0 0 1" armature="{p.J_rotor}"/>
        <geom class="vis" type="mesh" mesh="input_carrier" rgba=".62 .64 .67 1" density="7850"/>
        <geom class="vis" type="mesh" mesh="input_shaft" rgba=".62 .64 .67 1" density="7850"/>
        {rollers("input_roller", w.inp)}
      </body>
      <body name="ring_seat">
        <joint name="ring_axial" type="slide" axis="0 0 1" stiffness="{p.k_ax}" damping="20"/>
        <geom class="vis" type="sphere" size="0.001" mass="0.001"/>
        <body name="ring">
          <joint name="ring" type="ball"/>
          <geom class="vis" type="mesh" mesh="ring" rgba=".93 .45 .12 1" density="7850"/>
        </body>
      </body>
      <body name="output" quat="{' '.join(map(str, oq))}">
        <joint name="output" type="hinge" axis="0 0 1"/>
        <geom class="vis" type="mesh" mesh="output_carrier" rgba=".2 .4 .7 1" density="2700"/>
        {rollers("output_roller", w.out)}
        <body name="arm" pos="0 0 {(C.housings(None)[2]['top'] + 18) * 1e-3}">
          <geom class="vis" type="box" size=".006 {p.arm_len / 2} .006" pos="0 {p.arm_len / 2} 0" rgba=".3 .3 .32 1" density="2700"/>
          <geom name="tip" class="vis" type="sphere" size=".02" pos="0 {p.arm_len} 0" mass="{p.arm_mass}" rgba=".8 .1 .1 1"/>
        </body>
      </body>
    </body>
  </worldbody>
  <actuator><motor name="target_deg" joint="input" gear="0" ctrllimited="true" ctrlrange="-90 90"/></actuator>
</mujoco>"""


class Wobble2Sim:
    def __init__(self, w=None, p=None):
        self.w = w or D.Wobble2()
        self.p = p or Params()
        self.m = mujoco.MjModel.from_xml_string(make_xml(self.w, self.p, mesh_files()))
        self.d = mujoco.MjData(self.m)
        self.pre = D.preload_mm(self.w)
        self.K = np.array([D.k_contact(self.w.inp)] * self.w.inp.rollers + [D.k_contact(self.w.out)] * self.w.out.rollers)
        j = lambda n: self.m.joint(n)  # noqa: E731
        self.qin, self.vin = j("input").qposadr[0], j("input").dofadr[0]
        self.qout, self.vout = j("output").qposadr[0], j("output").dofadr[0]
        self.qax, self.vax = j("ring_axial").qposadr[0], j("ring_axial").dofadr[0]
        self.qr, self.vr = j("ring").qposadr[0], j("ring").dofadr[0]
        self.b_ring, self.b_in, self.b_out = (self.m.body(n).id for n in ("ring", "input", "output"))
        self.mode, self.target, self.volts, self.phi0 = "hold", 0.0, 0.0, 0.0
        self.i_motor = 0.0
        self.pin_el = C.pin_el(json.loads((OUT / "cad" / "parts.json").read_text())["ring_band_deg"]) \
            if (OUT / "cad" / "parts.json").exists() else 0.0
        self.m.dof_damping[self.vout] = self.p.out_damping
        self.reset(0.0, 0.0)

    def reset(self, phi, dphi):
        m, d, w = self.m, self.d, self.w
        mujoco.mj_resetData(m, d)
        d.qpos[self.qin], d.qvel[self.vin] = phi, dphi
        d.qpos[self.qout] = D.out_angle(w, phi) - D.out_angle(w, 0.0)
        d.qvel[self.vout] = dphi / w.ratio
        d.qpos[self.qr:self.qr + 4] = quat(D.ring_rot(w, phi))
        # ball joint velocity is in the ring's local frame
        d.qvel[self.vr:self.vr + 3] = D.ring_rot(w, phi).T @ (ring_rate(w, phi) * dphi)
        mujoco.mj_forward(m, d)

    def apply(self):
        m, d, w, p = self.m, self.d, self.w, self.p
        phi, dphi = d.qpos[self.qin], d.qvel[self.vin]
        c = D.contacts(w, phi)
        n, Cc, t = c["n"], c["C"] * 1e-3, c["t"]
        out = c["row"] == 1
        # deviations from the conjugate poses (CAD frame, rotations about the centre)
        Rring = G.T @ d.xmat[self.b_ring].reshape(3, 3)
        Rid = D.ring_rot(w, phi)
        r = rotvec(Rring @ Rid.T)
        s_ax = d.qpos[self.qax]
        om_ring = G.T @ self._omega(self.b_ring)
        om_id = ring_rate(w, phi) * dphi
        dout = d.qpos[self.qout] - (D.out_angle(w, phi) - D.out_angle(w, 0.0))
        dout_v = d.qvel[self.vout] - dphi / w.ratio
        # roller minus ring displacement at C (ring translates only along z on the slide)
        dr = -(np.array([0, 0, s_ax]) + np.cross(r, Cc))
        vr = -(np.array([0, 0, d.qvel[self.vax]]) + np.cross(om_ring - om_id, Cc))
        dr[out] += np.cross([0, 0, dout], Cc[out])
        vr[out] += np.cross([0, 0, dout_v], Cc[out])
        delta = self.pre - (n * dr).sum(1) * 1e3
        ddot = -(n * vr).sum(1)
        Q = np.where(delta > 0, self.K * np.maximum(delta, 0) ** (10 / 9) + p.c_contact * ddot, 0.0)
        Q = np.maximum(Q, 0)
        F = Q[:, None] * n - (w.mu_roll * Q * np.tanh(dphi / 0.5))[:, None] * t     # on each roller
        # ring: -F at C; pins on the bisector, with sliding friction
        M_ring = -np.cross(Cc, F).sum(0)
        F_ring = -F.sum(0)
        a = Rid @ np.array([0, 0, 1.0])
        b = (a + [0, 0, 1.0]) / np.linalg.norm(a + [0, 0, 1.0])
        M_pin = -(p.k_pin * (r @ b) + p.c_pin * ((om_ring - om_id) @ b)) * b
        Rp = (C.RHO_SEAT + 2) * 1e-3
        F_pin = np.linalg.norm(M_pin) / (2 * Rp)
        for az in (0.0, np.pi):
            P = Rring @ (Rp * np.array([np.cos(self.pin_el) * np.cos(az), np.cos(self.pin_el) * np.sin(az), np.sin(self.pin_el)]))
            v = np.cross(om_ring, P)
            sp = np.linalg.norm(v)
            if sp > 1e-9:
                f = -w.mu_pin * F_pin * np.tanh(sp / 0.002) * v / sp
                M_pin = M_pin + np.cross(P, f)
        self.ground_moment = -M_pin                                 # the pins push back on the housing
        self._wrench(self.b_ring, F_ring, M_ring + M_pin)
        # carriers: the rollers' forces, as torques about the axis (the hinges take the rest)
        inn = ~out
        tau_in = np.cross(Cc[inn], F[inn]).sum(0)[2]
        tau_out = np.cross(Cc[out], F[out]).sum(0)[2]
        # motor
        om = dphi
        if self.mode == "lock":
            i_cmd = -(p.k_lock * (phi - self.phi0) + p.c_lock * om) / p.Kt
        elif self.mode == "hold":
            ang = d.qpos[self.qout]
            T_cmd = -(p.kp * (ang - self.target) + p.kd * d.qvel[self.vout]) / w.ratio
            i_cmd = T_cmd / p.Kt
        elif self.mode == "speed":
            err = self.speed - om
            self.i_int = getattr(self, "i_int", 0.0) + err * p.dt
            i_cmd = 2.0 * err + 300.0 * self.i_int
        else:
            i_cmd = (self.volts - p.Kt * om) / p.R
        i_cmd = np.clip(i_cmd, (-p.V - p.Kt * om) / p.R, (p.V - p.Kt * om) / p.R)
        self.i_motor = float(np.clip(i_cmd, -p.I_max, p.I_max))
        tau_m = p.Kt * self.i_motor
        d.qfrc_applied[self.vin] = tau_m + tau_in
        d.qfrc_applied[self.vout] = tau_out - p.load_torque * np.tanh(d.qvel[self.vout] / 0.05)
        return dict(Q=Q, delta=delta, phi=phi, dphi=dphi, tau_m=tau_m, tau_in=tau_in, row=c["row"], k=c["k"])

    def _omega(self, b):
        v6 = np.zeros(6)
        mujoco.mj_objectVelocity(self.m, self.d, mujoco.mjtObj.mjOBJ_BODY, b, v6, 0)
        return v6[:3]

    def _wrench(self, b, F, M_about_origin):
        com = G.T @ self.d.xipos[b]
        orig = G.T @ self.d.xpos[self.m.body("frame").id]
        Mc = M_about_origin - np.cross(com - orig, F)
        self.d.xfrc_applied[b, :3] = G @ F
        self.d.xfrc_applied[b, 3:] = G @ Mc

    def step(self, n=1):
        for _ in range(n):
            mujoco.mj_step1(self.m, self.d)
            info = self.apply()
            mujoco.mj_step2(self.m, self.d)
        return info


# ---------------------------------------------------------------- scenarios

def arm_torque(sim):
    """gravity torque of the output's load about the CAD axis (N m)"""
    m, d = sim.m, sim.d
    b = sim.b_out
    mass = m.body_subtreemass[b]
    c = G.T @ d.subtree_com[b]
    F = G.T @ (mass * np.array([0, 0, -9.81]))
    return np.cross(c, F)[2]


def hold(w=None, p=None, T_end=0.15, mode="hold", every=100):
    sim = Wobble2Sim(w, p)
    sim.mode = mode
    tr = dict(t=[], ang=[], Qmax=[], i=[])
    info = None
    for k in range(int(T_end / sim.p.dt)):
        info = sim.step()
        if k % every == 0:
            tr["t"].append(sim.d.time)
            tr["ang"].append(np.degrees(sim.d.qpos[sim.qout]))
            tr["Qmax"].append(info["Q"].max())
            tr["i"].append(sim.i_motor)
        if not np.isfinite(sim.d.qpos).all():
            break
    return sim, info, {k: np.array(v) for k, v in tr.items()}


def selftest():
    p = Params(out_damping=40.0)
    sim, info, tr = hold(p=p, T_end=0.25, mode="lock")
    assert np.isfinite(sim.d.qpos).all(), "blew up"
    T = arm_torque(sim)
    s = D.solve(sim.w, info["phi"], T=T * 1e3)
    err = np.abs(info["Q"] - s["Q"]).max()
    print(f"model: masses (kg) ring {sim.m.body_mass[sim.b_ring]:.3f}, input {sim.m.body_subtreemass[sim.b_in]:.3f}, "
          f"output+arm {sim.m.body_subtreemass[sim.b_out]:.3f}")
    print(f"held {abs(T):.1f} N m: motor {sim.i_motor:.2f} A ({abs(sim.i_motor * p.Kt):.3f} N m, expected "
          f"{abs(T / sim.w.ratio):.3f}); roller loads {info['Q'].min():.0f}-{info['Q'].max():.0f} N")
    print(f"statics vs design.solve: worst roller-load difference {err:.1f} N (of {s['Q'].max():.0f} N)")
    assert err < 0.05 * s["Q"].max() + 10, "sim statics disagree with design.solve"
    print("selftest ok")


def ratio_run(dphi=150.0, T_end=0.06):
    p = replace(Params(), arm_mass=0.01)
    sim = Wobble2Sim(p=p)
    sim.reset(0.0, dphi)
    sim.mode, sim.volts = "volts", p.Kt * dphi
    ph, out = [], []
    for k in range(int(T_end / p.dt)):
        info = sim.step()
        if k % 50 == 0:
            ph.append(info["phi"])
            out.append(sim.d.qpos[sim.qout])
    ph, out = np.array(ph), np.array(out)
    return dict(phi=ph, out=out, ratio=1 / np.polyfit(ph, out, 1)[0], expected=sim.w.ratio)


def eff_run(speed=150.0, T_load=10.0, turns=2, w=None):
    """constant input speed (speed loop) against a constant dyno torque on the output: efficiency
    = output power / motor shaft power, averaged over whole input turns after settling"""
    p = replace(Params(), arm_mass=0.01, load_torque=T_load)
    sim = Wobble2Sim(w, p)
    sim.reset(0.0, speed)
    sim.mode, sim.speed = "speed", speed
    n0 = max(int(0.02 / p.dt), int(2 * np.pi / speed / p.dt))
    n1 = n0 + turns * int(2 * np.pi / speed / p.dt)
    Pin, Pout = [], []
    for k in range(n1):
        info = sim.step()
        if k >= n0:
            Pin.append(info["tau_m"] * info["dphi"])
            Pout.append(T_load * abs(sim.d.qvel[sim.vout]))
    return dict(eta=float(np.mean(Pout) / np.mean(Pin)), T=T_load, speed=speed, design=D.efficiency(sim.w, T_load))


def shake_run(dphi=314.0, T_end=0.03):
    p = replace(Params(), arm_mass=0.01)
    sim = Wobble2Sim(p=p)
    sim.reset(0.0, dphi)
    sim.mode, sim.volts = "volts", p.Kt * dphi
    t, H = [], []
    for k in range(int(T_end / p.dt)):
        sim.step()
        if k % 10 == 0:
            mujoco.mj_subtreeVel(sim.m, sim.d)
            t.append(sim.d.time)
            H.append(G.T @ sim.d.subtree_angmom[0].copy())
    t, H = np.array(t), np.array(H)
    dH = np.gradient(H, t, axis=0)
    M = np.hypot(dH[:, 0], dH[:, 1])
    s = t > 0.005
    return dict(t=t, M=M, rms=float(np.sqrt(np.mean(M[s] ** 2))))


def gif(path, revs=2, frames_per_rev=45, dphi=60.0):
    from PIL import Image, ImageDraw, ImageFont
    global CUTAWAY
    CUTAWAY = True
    p = replace(Params(), arm_mass=0.01, arm_len=0.12)
    sim = Wobble2Sim(p=p)
    CUTAWAY = False
    sim.reset(0.0, dphi)
    sim.mode, sim.volts = "volts", p.Kt * dphi
    r = mujoco.Renderer(sim.m, 540, 960)
    cam = mujoco.MjvCamera()
    cam.lookat[:] = (0.0, 0.0, 0.0)
    cam.distance, cam.azimuth, cam.elevation = 0.36, 35.0, -20.0
    per = int(2 * np.pi / frames_per_rev / dphi / p.dt)
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    ims = []
    for f in range(revs * frames_per_rev):
        sim.step(per)
        r.update_scene(sim.d, cam)
        im = Image.fromarray(r.render())
        dr = ImageDraw.Draw(im)
        dr.text((14, 10), f"wobble2 {abs(sim.w.ratio):.0f}:1 - one wobble ring, MuJoCo 4.0 roller contacts, slow motion",
                fill=(20, 20, 25), font=font)
        dr.text((14, 36), f"input {sim.d.qpos[sim.qin] / 2 / np.pi:5.2f} rev   output {np.degrees(sim.d.qpos[sim.qout]):7.2f} deg",
                fill=(60, 60, 70), font=font)
        ims.append(im)
    ims[0].save(path, save_all=True, append_images=ims[1:], duration=45, loop=0)


def demo():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    res = {}
    print("ratio ...")
    rt = ratio_run()
    res["ratio"] = dict(measured=rt["ratio"], expected=rt["expected"])
    print(f"  measured {rt['ratio']:.3f}, expected {rt['expected']:.3f}")
    print("hold ...")
    sim, info, tr = hold(p=Params(), T_end=0.25)
    T = arm_torque(sim)
    s = D.solve(sim.w, info["phi"], T=T * 1e3)
    res["hold"] = dict(T_Nm=float(abs(T)), Qmax=float(info["Q"].max()), solve_Qmax=float(s["Q"].max()), motor_A=sim.i_motor)
    print(f"  {abs(T):.1f} N m: roller max {info['Q'].max():.0f} N (design.solve {s['Q'].max():.0f}), motor {sim.i_motor:.2f} A")
    print("efficiency ...")
    ef = [eff_run(sp) for sp in (20.0, 60.0, 150.0, 314.0)]
    res["efficiency"] = ef
    for e in ef:
        print(f"  {e['T']:.0f} N m out at {e['speed'] * 60 / 2 / np.pi:.0f} rpm in: sim {e['eta'] * 100:.1f} %, "
              f"design.efficiency (quasi-static) {e['design'] * 100:.1f} %")
    print("shake at 3000 rpm input ...")
    sh = shake_run()
    res["shake"] = dict(rms_Nm=sh["rms"], wobble1_twin_rms=0.004, wobble1_single_rms=1.03)
    print(f"  frame moment rms {sh['rms']:.2f} N m (wobble1: twin 0.004, one block 1.03)")
    fig, ax = plt.subplots(1, 4, figsize=(19, 4))
    ax[3].plot([e["speed"] * 60 / 2 / np.pi for e in ef], [e["eta"] * 100 for e in ef], "o-", label="sim")
    ax[3].axhline(ef[0]["design"] * 100, ls="--", c="k", lw=0.8, label="quasi-static model")
    ax[3].set_xlabel("input speed (rpm)")
    ax[3].set_ylabel("efficiency (%)")
    ax[3].set_title(f"efficiency at {ef[0]['T']:.0f} N m out (dyno load)")
    ax[3].legend(fontsize=8)
    ax[0].plot(rt["phi"] / 2 / np.pi, np.degrees(rt["out"]), label="output")
    ax[0].plot(rt["phi"] / 2 / np.pi, np.degrees(rt["phi"] / rt["expected"]), "--", label=f"input / {rt['expected']:.0f}")
    ax[0].set_xlabel("input (rev)")
    ax[0].set_ylabel("output (deg)")
    ax[0].set_title(f"ratio: measured {rt['ratio']:.3f}")
    ax[0].legend()
    ax[1].plot(tr["t"] * 1e3, tr["Qmax"], label="max roller load (N)")
    ax[1].axhline(s["Q"].max(), ls="--", c="k", lw=0.8, label="design.solve")
    ax[1].set_xlabel("ms")
    ax[1].set_title(f"hold {abs(T):.1f} N m with the output servo")
    ax[1].legend(fontsize=8)
    ax[2].plot(sh["t"] * 1e3, sh["M"])
    ax[2].set_xlabel("ms")
    ax[2].set_ylabel("frame moment |Mxy| (N m)")
    ax[2].set_title(f"shaking at 3000 rpm in ({sh['rms']:.2f} N m rms)")
    fig.tight_layout()
    fig.savefig(OUT / "sim_demo.png", dpi=110)
    (OUT / "sim.json").write_text(json.dumps(res, indent=1, default=float))
    print("gif ...")
    gif(OUT / "sim.gif")
    print("wrote out/sim_demo.png, out/sim.json, out/sim.gif")


def run_interactive():
    import time
    import mujoco.viewer
    sim = Wobble2Sim()
    with mujoco.viewer.launch_passive(sim.m, sim.d) as v:
        while v.is_running():
            t0 = time.time()
            sim.target = np.radians(sim.d.ctrl[0])
            sim.step(200)
            v.sync()
            dt = time.time() - t0
            if dt < 1 / 60:
                time.sleep(1 / 60 - dt)


if __name__ == "__main__":
    {"selftest": selftest, "interactive": run_interactive, "demo": demo}[_MODE]()
