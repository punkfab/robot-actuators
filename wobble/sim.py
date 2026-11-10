"""
MuJoCo (4.0) model of the twin wobble reducer, from the CAD meshes and the design module.

The blocks and the output ring are FREE bodies: no joints, no bearing on the output. What holds
them is the roller contact law, evaluated every step from the same groove geometry the CAD cut:
for each of the 44 rollers, the approach of roller into wall (preload + the bodies' deviation
from their conjugate poses at the current input angle) gives a Palmgren/Hertz line-contact
force plus damping, along the wall normal, and rolling friction along the path. MuJoCo can't
collide a roller with a concave wavy groove (it would see the groove's convex hull), so the
contact is a custom force law applied through xfrc_applied (mujoco-sim skill: custom forces).
The input shaft is a hinge driven by a current- and voltage-limited motor; each wobble bearing
holds its block's tilt direction (as in design.solve).

    python wobble/sim.py              interactive viewer (slow motion: dt is 5 us)
    python wobble/sim.py --demo       scenarios -> out/sim_*.png, out/sim.gif
    python wobble/sim.py --selftest   build, hold a load, check against design.solve

Scenarios (--demo):
  ratio    spin the input; output/input speed must be 1/65; input torque ripple
  hold     an arm on the output, held by the motor; output deflection and roller loads compared
           with design.solve's statics for the same load (the known-physics check)
  ab       the same hold with ONE block and no output bearing: the output's free wobble lets the
           arm fall (design.py's finding), against the twin
  eff      lifting at speed: efficiency with rolling friction, and 100% without (check)
  shake    3000 rpm input, no load: frame reaction moment, one block (with its aux bearing)
           against the twin, whose opposed nutations cancel
"""

import os
import sys

_MODE = "demo" if "--demo" in sys.argv else "selftest" if "--selftest" in sys.argv else "interactive"
os.environ.setdefault("MUJOCO_GL", "glfw" if _MODE == "interactive" else "egl")

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

# world <- CAD: CAD x -> world z (up), CAD y -> world x, CAD z -> world y (the actuator axis lies
# horizontal, the output lug points along world x, gravity is CAD -x)
G = np.array([[0, 1, 0], [0, 0, 1], [1, 0, 0.0]])
CUTAWAY = False          # build the model with the -y half of the housings and output ring cut away (the gif)


@dataclass
class Params:
    dt: float = 5e-6                 # s   the roller contacts are stiff (~1.7e8 N/m each)
    mu_roll: float = 0.02            # rolling friction coefficient per roller (stem in a plain bore, assumed)
    c_contact: float = 300.0         # N s/m damping per roller contact
    k_in: float = 2e6                # N m/rad wobble-bearing tilt stiffness (holds the block's tilt direction)
    c_in: float = 20.0               # N m s/rad
    k_aux: float = 0.0               # N/m output aux radial bearing (the single-block variant needs it)
    h_aux: float = 0.025             # m   above the block centre
    k_aux_tilt: float = 0.0          # N m/rad aux bearing tilt stiffness (a crossed roller main bearing)
    c_aux_tilt: float = 50.0         # N m s/rad
    arm_mass: float = 4.0            # kg  load at the arm tip
    arm_len: float = 0.30            # m   tip distance from the axis
    arm_offset: float = 0.0          # m   tip offset along the axis (a tilting moment on the output)
    # motor: a stand-in 48 V pancake BLDC (assumed numbers)
    Kt: float = 0.12                 # N m/A
    R: float = 0.20                  # ohm
    V: float = 48.0                  # V
    I_max: float = 15.0              # A    -> 1.8 N m
    J_rotor: float = 2e-4            # kg m^2
    kp: float = 3000.0               # N m/rad on the output angle (motor current command, output side)
    kd: float = 120.0                # N m s/rad
    k_lock: float = 50.0             # N m/rad input held in place ("lock" mode, the statics check)
    out_damping: float = 0.0         # N m s/rad on the output's free joint (implicit; settles the statics check)


def rotvec(Rm):
    """small-rotation vector of a rotation matrix (exact via the axis-angle formula)"""
    c = np.clip((np.trace(Rm) - 1) / 2, -1, 1)
    a = np.arccos(c)
    v = np.array([Rm[2, 1] - Rm[1, 2], Rm[0, 2] - Rm[2, 0], Rm[1, 0] - Rm[0, 1]])
    return v / 2 if a < 1e-9 else v * a / (2 * np.sin(a))


def quat(Rm):
    q = np.zeros(4)
    mujoco.mju_mat2Quat(q, Rm.flatten())
    return q


# ---------------------------------------------------------------- model

def meshes(w):
    """write the meshes the model uses (CAD frame, mm) into out/sim/"""
    SIM.mkdir(parents=True, exist_ok=True)
    cad = OUT / "cad"
    if not (cad / "output_ring.stl").exists():
        C.main()
    names = {}
    for n in ("housing_a", "housing_b", "output_ring", "input_shaft", "clevis", "motor", "block_canonical",
              "roller_canonical", "wobble_bearing"):
        names[n] = str(cad / f"{n}.stl")
    if CUTAWAY:
        box = trimesh.creation.box(extents=(400, 400, 400))
        box.apply_translation((0, -200, 0))
        for n in ("housing_a", "housing_b", "output_ring"):
            cut = C.boolean(trimesh.load(names[n]), box, "diff")
            cut.export(SIM / f"{n}_cut.stl")
            names[n] = str(SIM / f"{n}_cut.stl")
        empty = trimesh.creation.box(extents=(0.1, 0.1, 0.1))
        empty.apply_translation((0, 0, -500))
        empty.export(SIM / "hidden.stl")
        names["clevis"] = str(SIM / "hidden.stl")
    # block 2 is the mirror image of block 1
    b = trimesh.load(names["block_canonical"])
    A, _ = D.placements(D.Wobble())[1] if len(w.blocks) > 1 else (np.eye(3), None)
    M = np.eye(4)
    M[:3, :3] = A
    b.apply_transform(M)
    b.export(SIM / "block_b.stl")
    names["block_b"] = str(SIM / "block_b.stl")
    return names


def make_xml(w, p, files):
    pl = D.placements(w)
    g = np.array([0, 0, -9.81])                       # world; CAD -x
    asset = "\n".join(f'<mesh name="{n}" file="{f}" scale="0.001 0.001 0.001" inertia="exact"/>' for n, f in files.items())
    blocks = []
    for k, (A, o) in enumerate(pl):
        Rb = A @ D.block_rot(w, 0.0) @ A.T
        rollers = []
        for row in (w.fixed, w.output):
            u, _ = D.axes(w, row)
            for uj in u:
                Rr = C.rot_to(A @ uj)
                rollers.append(f'<geom class="vis" type="mesh" mesh="roller_canonical" quat="{" ".join(map(str, quat(Rr)))}" '
                               f'rgba=".75 .77 .8 1" density="7850"/>')
        mesh = "block_canonical" if k == 0 else "block_b"
        blocks.append(f'''
      <body name="block{k}" pos="{' '.join(map(str, G @ o * 1e-3))}" quat="{' '.join(map(str, quat(G @ Rb)))}">
        <freejoint name="block{k}"/>
        <geom class="vis" type="mesh" mesh="{mesh}" rgba=".93 .45 .12 1" density="7850"/>
        <geom class="vis" type="mesh" mesh="wobble_bearing" rgba=".85 .86 .9 1" density="0.001"/>
        {''.join(rollers)}
      </body>''')
    arm_len = p.arm_len
    return f"""<mujoco model="wobble">
  <compiler angle="radian"/>
  <option timestep="{p.dt}" integrator="Euler" gravity="{g[0]} {g[1]} {g[2]}"/>
  <visual><global offwidth="1280" offheight="720"/><headlight ambient=".4 .4 .4"/></visual>
  <default><default class="vis"><geom contype="0" conaffinity="0"/></default></default>
  <asset>{asset}
    <texture name="sky" type="skybox" builtin="gradient" rgb1=".97 .97 .98" rgb2=".78 .79 .82" width="64" height="64"/>
  </asset>
  <worldbody>
    <light pos="0 -0.4 0.6" dir="0 .5 -1" diffuse=".8 .8 .8"/>
    <body name="frame" quat="{' '.join(map(str, quat(G)))}">
      <geom class="vis" type="mesh" mesh="housing_a" rgba=".8 .81 .83 .35"/>
      <geom class="vis" type="mesh" mesh="housing_b" rgba=".8 .81 .83 .35"/>
      <geom class="vis" type="mesh" mesh="clevis" rgba=".15 .15 .16 .25"/>
      <geom class="vis" type="mesh" mesh="motor" rgba=".12 .12 .13 1"/>
    </body>
    <body name="input" quat="{' '.join(map(str, quat(G)))}">
      <joint name="input" type="hinge" axis="0 0 1" armature="{p.J_rotor}"/>
      <geom class="vis" type="mesh" mesh="input_shaft" rgba=".6 .62 .65 1" density="7850"/>
    </body>
    {''.join(blocks)}
    <body name="output" quat="{' '.join(map(str, quat(G)))}">
      <freejoint name="output"/>
      <geom class="vis" type="mesh" mesh="output_ring" rgba=".2 .4 .7 .5" density="7850"/>
      <body name="arm" pos="0 {(C.RING_R + C.LUG[0]) * 1e-3} 0">
        <geom class="vis" type="box" size=".01 {arm_len / 2} .006" pos="0 {arm_len / 2} 0" rgba=".3 .3 .32 1" density="2700"/>
        <geom class="vis" type="box" size=".006 .006 {abs(p.arm_offset) / 2 + 1e-4}" pos="0 {arm_len} {p.arm_offset / 2}" rgba=".3 .3 .32 1" density="2700"/>
        <geom name="tip" class="vis" type="sphere" size=".02" pos="0 {arm_len} {p.arm_offset}" mass="{p.arm_mass}" rgba=".8 .1 .1 1"/>
      </body>
    </body>
  </worldbody>
  <actuator>
    <motor name="target_deg" joint="input" gear="0" ctrllimited="true" ctrlrange="-90 90"/>
  </actuator>
</mujoco>"""


class Wobble:
    """model + per-step roller contact forces"""

    def __init__(self, w=None, p=None, xml_only=False):
        self.w = w or D.Wobble()
        self.p = p or Params()
        files = meshes(self.w)
        self.m = mujoco.MjModel.from_xml_string(make_xml(self.w, self.p, files))
        self.d = mujoco.MjData(self.m)
        self.pl = D.placements(self.w)
        self.nb = len(self.pl)
        self.K = D.k_contact(self.w)                  # N / mm^(10/9)
        self.pre = self.w.preload                      # mm
        self.bid = [self.m.body(f"block{k}").id for k in range(self.nb)] + [self.m.body("output").id]
        self.jin = self.m.joint("input")
        self.qin, self.vin = self.jin.qposadr[0], self.jin.dofadr[0]
        self.i_motor = 0.0
        self.log = {}
        self.ground = np.zeros(6)                      # reaction force/moment on the frame (CAD frame, about origin)
        self.target = 0.0
        self.mode = "hold"                             # "hold" (PD on the output angle) or "volts"
        self.volts = 0.0
        self.phi0 = 0.0
        self.reset(0.0, 0.0)
        if self.p.out_damping:
            va = self.m.jnt_dofadr[self.m.body_jntadr[self.bid[-1]]]
            self.m.dof_damping[va:va + 6] = self.p.out_damping

    # ---- kinematics of the ideal (conjugate) motion, CAD frame
    def ideal(self, phi, dphi):
        w = self.w
        Rb = D.block_rot(w, phi)
        Om = D.block_spin_axis(w, phi)
        out = []
        for A, o in self.pl:
            out.append((A @ Rb @ A.T, o * 1e-3, np.linalg.det(A) * (A @ Om) * dphi))
        out.append((D.Rz(D.output_angle(w, phi)), np.zeros(3), np.array([0, 0, dphi / w.ratio])))
        return out

    def reset(self, phi, dphi):
        m, d = self.m, self.d
        mujoco.mj_resetData(m, d)
        d.qpos[self.qin], d.qvel[self.vin] = phi, dphi
        for b, (Rs, ps, ws) in zip(self.bid, self.ideal(phi, dphi)):
            j = m.body_jntadr[b]
            qa, va = m.jnt_qposadr[j], m.jnt_dofadr[j]
            d.qpos[qa:qa + 3] = G @ ps
            d.qpos[qa + 3:qa + 7] = quat(G @ Rs)
            d.qvel[va:va + 3] = np.cross(G @ ws, G @ ps - G @ ps)      # body origin on the ideal path: no translation
            d.qvel[va + 3:va + 6] = (G @ Rs).T @ (G @ ws)               # free joint angular velocity is local
        mujoco.mj_forward(m, d)

    def body_state(self, b):
        """CAD-frame pose and velocity of a body's frame origin"""
        d = self.d
        R = G.T @ d.xmat[b].reshape(3, 3)
        x = G.T @ d.xpos[b]
        v6 = np.zeros(6)
        mujoco.mj_objectVelocity(self.m, d, mujoco.mjtObj.mjOBJ_BODY, b, v6, 0)
        return R, x, G.T @ v6[3:], G.T @ v6[:3]          # R, origin, v_origin, omega

    def output_angle(self):
        R = G.T @ self.d.xmat[self.bid[-1]].reshape(3, 3)
        return np.arctan2(R[1, 0], R[0, 0])

    # ---- forces
    def apply(self):
        m, d, w, p = self.m, self.d, self.w, self.p
        phi, dphi = d.qpos[self.qin], d.qvel[self.vin]
        ca = D.contact_arrays(w, phi)
        Cm = ca["C"] * 1e-3                                              # m
        ide = self.ideal(phi, dphi)
        st = [self.body_state(b) for b in self.bid]
        # deviation field d(X) = v + r x X and its rate, per body
        dev, rate = [], []
        for (R, x, vo, om), (Rs, ps, ws) in zip(st, ide):
            r = rotvec(R @ Rs.T)
            dev.append((x - ps - np.cross(r, ps), r))
            rate.append((vo, om, x, ws, ps))

        def disp(i, X):
            v, r = dev[i]
            return v + np.cross(r, X)

        def vel(i, X):
            vo, om, x, ws, ps = rate[i]
            return vo + np.cross(om, X - x) - np.cross(ws, X - ps)

        nb = self.nb
        blk, out = ca["block"], ca["out"]
        drel = np.zeros_like(Cm)
        vrel = np.zeros_like(Cm)
        for k in range(nb):
            s = blk == k
            drel[s] = disp(k, Cm[s])
            vrel[s] = vel(k, Cm[s])
        so = out
        drel[so] -= disp(nb, Cm[so])
        vrel[so] -= vel(nb, Cm[so])
        n, t = ca["n"], ca["t"]
        delta = self.pre - (n * drel).sum(1) * 1e3                      # mm
        ddot = -(n * vrel).sum(1)                                        # m/s
        Q = np.where(delta > 0, self.K * np.maximum(delta, 0) ** (10 / 9) + p.c_contact * ddot, 0.0)
        Q = np.maximum(Q, 0)
        fr = -p.mu_roll * Q * np.tanh(dphi / 0.5)
        F = Q[:, None] * n + fr[:, None] * t                              # on the roller (block)
        # wrenches about each body's COM, then to world
        wrench = [np.zeros(6) for _ in self.bid]
        com = [G.T @ d.xipos[b] for b in self.bid]

        def add(i, Fs, Xs):
            wrench[i][:3] += Fs.sum(0)
            wrench[i][3:] += np.cross(Xs - com[i], Fs).sum(0)

        for k in range(nb):
            s = blk == k
            add(k, F[s], Cm[s])
        add(nb, -F[so], Cm[so])
        fixed = ~so
        self.ground[:3] = -F[fixed].sum(0)
        self.ground[3:] = -np.cross(Cm[fixed], F[fixed]).sum(0)
        # wobble bearings: hold each block's tilt direction; the reaction is the input torque
        tau_in = 0.0
        Om = D.block_spin_axis(w, phi)
        for k, ((Rs, ps, ws), (R, x, vo, om)) in enumerate(zip(ide[:nb], st[:nb])):
            A = self.pl[k][0]
            ab = Rs @ np.array([0, 0, 1.0])
            P = np.eye(3) - np.outer(ab, ab)
            r = dev[k][1]
            Min = -p.k_in * P @ r - p.c_in * P @ (om - ws)
            wrench[k][3:] += Min
            self.ground[3:] -= Min                                        # through the shaft's bearings
            tau_in += Min @ (np.linalg.det(A) * A @ Om)
        # the single-block variant's aux radial bearing on the output
        if p.k_aux:
            R, x, vo, om = st[nb]
            X = np.array([0, 0, p.h_aux])
            u = disp(nb, X)
            Fa = -p.k_aux * np.array([u[0], u[1], 0])
            add(nb, Fa[None], X[None])
            self.ground[:3] -= Fa
            self.ground[3:] -= np.cross(X, Fa)
        if p.k_aux_tilt:
            r = dev[nb][1]
            om = st[nb][3] - ide[nb][2]
            Mt = -p.k_aux_tilt * np.array([r[0], r[1], 0]) - p.c_aux_tilt * np.array([om[0], om[1], 0])
            wrench[nb][3:] += Mt
            self.ground[3:] -= Mt
        for i, b in enumerate(self.bid):
            d.xfrc_applied[b, :3] = G @ wrench[i][:3]
            d.xfrc_applied[b, 3:] = G @ wrench[i][3:]
        # motor
        om_in = dphi
        if self.mode == "lock":
            i_cmd = -(p.k_lock * (phi - self.phi0) + 0.05 * om_in) / p.Kt
        elif self.mode == "hold":
            T_cmd = -(p.kp * (self.output_angle() - self.target) + p.kd * (om_in / w.ratio)) / w.ratio
            i_cmd = T_cmd / p.Kt
        else:
            i_cmd = (self.volts - p.Kt * om_in) / p.R
        i_hi = (p.V - p.Kt * om_in) / p.R
        i_lo = (-p.V - p.Kt * om_in) / p.R
        self.i_motor = float(np.clip(np.clip(i_cmd, i_lo, i_hi), -p.I_max, p.I_max))
        tau_m = p.Kt * self.i_motor
        d.qfrc_applied[self.vin] = tau_m - tau_in
        return dict(Q=Q, delta=delta, tau_in=tau_in, tau_m=tau_m, phi=phi, dphi=dphi, k=ca["k"], out=out, blk=blk)

    def step(self, n=1):
        """mj_step1 (kinematics for the current state), our forces, mj_step2: forces computed from
        the state they act on. (After a plain mj_step, xpos/xmat still describe the previous state,
        and a one-step-stale stiff force goes unstable well inside the stability limit.)"""
        for _ in range(n):
            mujoco.mj_step1(self.m, self.d)
            info = self.apply()
            mujoco.mj_step2(self.m, self.d)
        return info


# ---------------------------------------------------------------- scenarios

def arm_load(sim):
    """static load on the output from the arm's weight (CAD frame, about the origin): T, Fr, M, Fa"""
    m, d = sim.m, sim.d
    b = m.body("arm").id
    mass = m.body_subtreemass[b]
    c = G.T @ d.subtree_com[b]
    F = G.T @ (mass * np.array([0, 0, -9.81]))
    Mo = np.cross(c, F)
    # the output ring's own weight too
    bo = sim.bid[-1]
    Fo = G.T @ (m.body_mass[bo] * np.array([0, 0, -9.81]))
    return F + Fo, Mo + np.cross(G.T @ d.xipos[bo], Fo)


def hold(w=None, p=None, T_end=0.25, every=200, mode="hold"):
    sim = Wobble(w, p)
    sim.mode = mode
    tr = dict(t=[], ang=[], tilt=[], i=[], Qmax=[])
    for k in range(int(T_end / sim.p.dt)):
        info = sim.step()
        if k % every == 0:
            R = G.T @ sim.d.xmat[sim.bid[-1]].reshape(3, 3)
            tr["t"].append(sim.d.time)
            tr["ang"].append(np.degrees(sim.output_angle()))
            tr["tilt"].append(np.degrees(np.arccos(np.clip(R[2, 2], -1, 1))) * 60)
            tr["i"].append(sim.i_motor)
            tr["Qmax"].append(info["Q"].max())
            if not np.isfinite(sim.d.qpos).all() or abs(tr["ang"][-1]) > 60:
                break
    return sim, info, {k: np.array(v) for k, v in tr.items()}


def statics_check(sim, info):
    """compare the settled roller loads with design.solve under the same load"""
    w = sim.w
    F, Mo = arm_load(sim)
    phi = info["phi"]
    s = solve_general(w, phi, F, Mo * 1e3, sim.p.k_in * 1e3)
    err = np.abs(info["Q"] - s["Q"]).max()
    return s, err


def solve_general(w, phi, F, M, k_in):
    """design.solve with an arbitrary force F (N, at the origin) and moment M (N mm) on the output"""
    cs = D.contacts(w, phi)
    nb = len(w.blocks)
    io = 6 * nb
    J = np.array([c["J"] for c in cs])
    K = D.k_contact(w)
    pl = D.placements(w)
    a_b = D.block_rot(w, phi) @ np.array([0, 0, 1.0])
    nd = 6 * (nb + 1)
    Pin = np.zeros((nd, nd))
    for bi, (A, o) in enumerate(pl):
        ab = A @ a_b
        Pin[6 * bi + 3:6 * bi + 6, 6 * bi + 3:6 * bi + 6] = np.eye(3) - np.outer(ab, ab)
    Gv = np.zeros(nd)
    Gv[io:io + 3] = F
    Gv[io + 3:] = M
    x = np.zeros(nd)
    for _ in range(100):
        d = np.maximum(w.preload + J @ x, 0)
        Q = K * d ** (10 / 9)
        g = J.T @ Q - Gv + k_in * Pin @ x
        H = (J * (K * 10 / 9 * d ** (1 / 9))[:, None]).T @ J + k_in * Pin + 1e-6 * np.eye(nd)
        dx = -np.linalg.solve(H, g)
        x += dx
        if np.abs(dx).max() < 1e-12:
            break
    d = np.maximum(w.preload + J @ x, 0)
    return dict(Q=K * d ** (10 / 9), x=x, xo=x[io:])


def selftest():
    p = Params(out_damping=60.0)
    sim, info, tr = hold(p=p, T_end=0.12, mode="lock")
    assert np.isfinite(sim.d.qpos).all(), "blew up"
    s, err = statics_check(sim, info)
    print(f"model: {sim.m.nbody} bodies; masses (kg): " +
          ", ".join(f"{sim.m.body(b).name} {sim.m.body_mass[b]:.3f}" for b in sim.bid) +
          f", arm+tip {sim.m.body_subtreemass[sim.m.body('arm').id]:.2f}")
    print(f"held: output {tr['ang'][-1]:+.3f} deg, tilt {tr['tilt'][-1]:.2f} arcmin, motor {sim.i_motor:.2f} A; "
          f"roller loads {info['Q'].min():.0f}-{info['Q'].max():.0f} N")
    print(f"statics vs design.solve: worst roller-load difference {err:.1f} N (of {s['Q'].max():.0f} N)")
    assert err < 0.05 * s["Q"].max() + 10, "sim statics disagree with design.solve"
    print("selftest ok")


def run_interactive():
    import mujoco.viewer
    sim = Wobble()
    sim.d.ctrl[0] = 0.0
    with mujoco.viewer.launch_passive(sim.m, sim.d) as v:
        import time
        while v.is_running():
            t0 = time.time()
            sim.target = np.radians(sim.d.ctrl[0])
            sim.step(400)                                  # 2 ms of sim per frame: slow motion
            v.sync()
            dt = time.time() - t0
            if dt < 1 / 60:
                time.sleep(1 / 60 - dt)


if __name__ == "__main__":
    if _MODE == "selftest":
        selftest()
    elif _MODE == "interactive":
        run_interactive()
    else:
        import demo
        demo.main()
