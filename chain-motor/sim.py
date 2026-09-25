"""
MuJoCo model of the 2-CELL VARIABLE CHAIN MOTOR — dynamic checks + a viewer.

The CAD meshes ride a real joint tree, and the GEARS are equality constraints (not contacts):

    world ─ thA (hinge)  rotor A (output)            motor A torque: rotor A vs world (= A stator)
    world ─ phi (hinge)  swing link
                └ sB (hinge)  stator B               gear:  sB = r_s·phi         (joint equality)
                     └ rB (hinge)  rotor B           gear:  rB + r_r·thA − (r_r−r_s)·phi = 0
                                                            (fixed-tendon equality)
                                                     motor B torque: rotor B vs stator B (on rB)

rB is B's rotor-to-stator angle, i.e. what its commutation sees — so the constraints are exactly
chain.py's kinematics. `check` then measures, with dynamics rather than algebra:

  1. hinge neutrality  stall (output + hinge held), motors torqued: torque on the hinge?
                       (r_s = r_r vs 33:34 — should be T_B·(r_r − r_s))
  2. torque summing    output accel with B's leads swapped (adds) vs not swapped (cancels)
  3. hinge inertia     output held by a servo while the hinge swings: holding torque vs the
                       predicted −2·I_rotorB·φ̈ (rotor B turns 2φ in the world — the one coupling
                       the static argument misses)

    ../.venv/bin/python chain-motor/sim.py            # viewer: hinge swings while the motor spins
    ../.venv/bin/python chain-motor/sim.py check      # the three dynamic checks
    ../.venv/bin/python chain-motor/sim.py render     # headless -> out/sim.gif + montage
"""

import sys
import time
from math import pi
from pathlib import Path

import numpy as np
import mujoco

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from chain import ChainParams
from cad import Cad

OUT = HERE / "out"
P = ChainParams()
C = Cad(P)
KT = P.cell.kt
I_ROTOR = 2.5e-6          # 2204 bell + magnets + rotor gear about the axis [kg·m²]
M_ROTOR = 0.016


def build_xml(r_s=None, hinge_damping=0.0, servo_output=False, hinge_servo=False):
    r_r = P.r_r
    r_s = P.r_s if r_s is None else r_s
    m = 1e-3
    a = C.a * m
    ph = C.tooth_phase_deg * pi / 180
    q_ph = f"{np.cos(ph/2):.6f} 0 0 {np.sin(ph/2):.6f}"
    marker = lambda z, rgba: (f'<geom type="box" size="0.0035 0.0012 0.0012" pos="0.012 0 {z*m:.4f}" '
                              f'rgba="{rgba}" contype="0" conaffinity="0" mass="0"/>')
    zr = C.z_rotor + C.gear_t + 0.6
    zs = C.z_stator + C.gear_t + 0.6
    coup = f'<joint joint="phi" coef="{-(r_r - r_s):.8f}"/>' if abs(r_r - r_s) > 0 else ""
    acts = ('<motor name="mA" joint="thA" gear="1"/>'
            '<motor name="mB" joint="rB" gear="1"/>')
    if servo_output:
        acts += '<position name="holdA" joint="thA" kp="5" kv="0.02"/>'
    if hinge_servo:
        acts += '<position name="swing" joint="phi" kp="200" kv="4"/>'
    xml = f"""<mujoco model="chain_motor_2cell">
  <compiler meshdir="out" angle="radian" autolimits="true"/>
  <option timestep="0.0001" gravity="0 0 0" integrator="implicitfast"/>
  <visual><global offwidth="960" offheight="720"/>
    <headlight ambient="0.45 0.45 0.45" diffuse="0.55 0.55 0.55" specular="0.1 0.1 0.1"/></visual>
  <default><geom contype="0" conaffinity="0" density="1240"/>
    <equality solref="0.0005 1"/></default>
  <asset>
    <mesh name="base" file="base.stl" scale="{m} {m} {m}"/>
    <mesh name="link_lower" file="link_lower.stl" scale="{m} {m} {m}"/>
    <mesh name="link_upper" file="link_upper.stl" scale="{m} {m} {m}"/>
    <mesh name="stator_gear_A" file="stator_gear_A.stl" scale="{m} {m} {m}"/>
    <mesh name="stator_gear_B" file="stator_gear_B.stl" scale="{m} {m} {m}"/>
    <mesh name="rotor_gear" file="rotor_gear.stl" scale="{m} {m} {m}"/>
    <mesh name="motor" file="motor_2204.stl" scale="{m} {m} {m}"/>
    <mesh name="dowel" file="dowel.stl" scale="{m} {m} {m}"/>
    <texture type="skybox" builtin="gradient" rgb1="0.93 0.94 0.96" rgb2="0.72 0.75 0.80" width="64" height="64"/>
    <material name="pla" rgba="0.80 0.80 0.84 1"/>
    <material name="link" rgba="0.35 0.55 0.85 0.9"/>
    <material name="stator" rgba="0.25 0.40 0.85 1"/>
    <material name="rotor" rgba="0.85 0.30 0.25 1"/>
    <material name="motor" rgba="0.25 0.25 0.28 1"/>
  </asset>
  <worldbody>
    <light pos="0.05 -0.1 0.3" dir="-0.1 0.3 -1"/>
    <body name="aim" pos="0.0 0 0.018"/>
    <camera name="iso" pos="0.03 -0.10 0.12" mode="targetbody" target="aim"/>
    <geom type="mesh" mesh="base" material="pla"/>
    <geom type="mesh" mesh="stator_gear_A" material="stator"/>
    <geom type="mesh" mesh="motor" material="motor"/>
    {marker(zs, "1 1 0.2 1")}
    <body name="rotorA">
      <joint name="thA" type="hinge" axis="0 0 1"/>
      <inertial pos="0 0 {C.z_rotor*m:.4f}" mass="{M_ROTOR}" diaginertia="{I_ROTOR/2} {I_ROTOR/2} {I_ROTOR}"/>
      <geom type="mesh" mesh="rotor_gear" material="rotor" mass="0"/>
      <geom type="mesh" mesh="dowel" material="motor" mass="0"/>
      {marker(zr, "1 1 1 1")}
    </body>
    <body name="link">
      <joint name="phi" type="hinge" axis="0 0 1" damping="{hinge_damping}"/>
      <geom type="mesh" mesh="link_lower" material="link"/>
      <geom type="mesh" mesh="link_upper" material="link"/>
      <body name="statorB" pos="{a:.5f} 0 0" quat="{q_ph}">
        <joint name="sB" type="hinge" axis="0 0 1"/>
        <geom type="mesh" mesh="stator_gear_B" material="stator"/>
        <geom type="mesh" mesh="motor" material="motor"/>
        {marker(zs, "1 1 0.2 1")}
        <body name="rotorB">
          <joint name="rB" type="hinge" axis="0 0 1"/>
          <inertial pos="0 0 {C.z_rotor*m:.4f}" mass="{M_ROTOR}" diaginertia="{I_ROTOR/2} {I_ROTOR/2} {I_ROTOR}"/>
          <geom type="mesh" mesh="rotor_gear" material="rotor" mass="0"/>
          <geom type="mesh" mesh="dowel" material="motor" mass="0"/>
      <geom type="mesh" mesh="dowel" material="motor" mass="0"/>
          {marker(zr, "1 1 1 1")}
        </body>
      </body>
    </body>
  </worldbody>
  <tendon><fixed name="rotor_mesh"><joint joint="rB" coef="1"/><joint joint="thA" coef="{r_r}"/>{coup}</fixed></tendon>
  <equality>
    <joint name="stator_mesh" joint1="sB" joint2="phi" polycoef="0 {r_s} 0 0 0"/>
    <tendon name="rotor_mesh" tendon1="rotor_mesh"/>
  </equality>
  <actuator>{acts}</actuator>
</mujoco>"""
    path = HERE / "viewer.xml"
    path.write_text(xml)
    return path


def _model(**kw):
    m = mujoco.MjModel.from_xml_path(str(build_xml(**kw)))
    return m, mujoco.MjData(m)


def _q(m, d, j):
    return d.qpos[m.jnt_qposadr[m.joint(j).id]]


def _v(m, d, j):
    return d.qvel[m.jnt_dofadr[m.joint(j).id]]


def check():
    i_cmd = 2.0                                 # A of synchronized phase current
    tau = KT * i_cmd
    print(f"\n=== 2-cell chain motor — MuJoCo dynamic checks (i = {i_cmd} A, τ/cell = {tau*1e3:.1f} mN·m) ===")
    ok = True

    print("\n[1] hinge neutrality (stall, like the paper's Kt test): output AND hinge held by servos,")
    print("    both motors torqued — torque the hinge servo must supply:")
    for r_s, lab in ((1.0, "stator 34:34 (= rotor)"), (33 / 34, "stator 33:34 (1 tooth off)")):
        m, d = _model(r_s=r_s, servo_output=True, hinge_servo=True)
        d.ctrl[:] = [tau, -tau, 0.0, 0.0]       # B's leads swapped: torque on rB is −τ
        while d.time < 0.3:
            mujoco.mj_step(m, d)
        hinge = d.actuator_force[3]
        pred = tau * (P.r_r - r_s)
        print(f"   {lab:28s} hinge torque {hinge*1e3:7.3f} mN·m   (virtual work T_B·(r_r−r_s) = "
              f"{pred*1e3:.3f})   output stall {-d.actuator_force[2]*1e3:5.1f} mN·m")
        if r_s == 1.0:
            ok &= abs(hinge) < 0.02 * tau

    print("\n[2] torque summing at the output (hinge held, 20 ms spin-up from rest):")
    I_eff = I_ROTOR * (1 + P.r_r ** 2)
    for sign, lab in ((-1, "B leads swapped"), (+1, "B leads NOT swapped")):
        m, d = _model(hinge_damping=1e3)
        d.ctrl[:] = [tau, sign * tau]
        while d.time < 0.02:
            mujoco.mj_step(m, d)
        alpha = _v(m, d, "thA") / d.time
        print(f"   {lab:22s} effective output torque {alpha*I_eff*1e3:6.2f} mN·m  "
              f"(ideal sum {2*tau*1e3:.2f}, cancel 0)")
        if sign < 0:
            ok &= abs(alpha * I_eff / (2 * tau) - 1) < 0.05

    print("\n[3] hinge inertia: output HELD by a servo, hinge swung 0→90° in 0.1 s (cosine ramp)")
    m, d = _model(servo_output=True, hinge_servo=True)
    T = 0.1
    samples = []
    while d.time < T:
        s = d.time / T
        d.ctrl[2] = 0.0                                       # hold output at 0
        d.ctrl[3] = (pi / 2) * 0.5 * (1 - np.cos(pi * s))     # swing command
        mujoco.mj_step(m, d)
        samples.append((d.time, _q(m, d, "phi"), d.actuator_force[2]))
    t, phi, hold = map(np.array, zip(*samples))
    phidd = np.gradient(np.gradient(phi, t), t)
    pred = -2 * I_ROTOR * P.r_r * phidd * P.r_s            # −2·I_B·φ̈  (θ held)
    k = slice(len(t) // 10, -len(t) // 10)
    corr = np.corrcoef(hold[k], pred[k])[0, 1]
    scale = np.dot(hold[k], pred[k]) / np.dot(pred[k], pred[k])
    print(f"   peak holding torque {np.abs(hold[k]).max()*1e3:.2f} mN·m   predicted −2·I_B·φ̈ peak "
          f"{np.abs(pred[k]).max()*1e3:.2f} mN·m   (fit ×{scale:.2f}, r = {corr:.3f})")
    print(f"   output drift under hold {np.degrees(np.abs([_q(m, d, 'thA')])).max():.2f}°")
    ok &= corr > 0.95 and abs(scale - 1) < 0.15
    print(f"\nVERDICT: {'VALID' if ok else 'CHECK'} — the hinge is torque-neutral statically, B must be "
          "phase-reversed, and fast hinge motion costs the output servo 2·I_B·φ̈.")
    return ok


def _pose(m, d, phi, theta):
    """Kinematic pose satisfying both gear constraints."""
    d.qpos[m.jnt_qposadr[m.joint("phi").id]] = phi
    d.qpos[m.jnt_qposadr[m.joint("thA").id]] = theta
    d.qpos[m.jnt_qposadr[m.joint("sB").id]] = P.r_s * phi
    d.qpos[m.jnt_qposadr[m.joint("rB").id]] = -P.r_r * theta + (P.r_r - P.r_s) * phi
    mujoco.mj_forward(m, d)


def _traj(s):
    """s in [0,1): hinge swings 0→155°→0 while the rotor turns steadily."""
    phi = np.radians(155) * 0.5 * (1 - np.cos(2 * pi * s))
    theta = 2 * pi * 1.5 * s
    return phi, theta


def view(period_s=8.0):
    from mujoco import viewer as mjv
    m, d = _model()
    dt = 1 / 60
    print("\nviewer: hinge swings 0→155° while the rotors spin — yellow = stator marks, white = rotor"
          "\n  marks; B's white mark keeps a fixed angle to its yellow one (commutation undisturbed)")
    with mjv.launch_passive(m, d) as v:
        t0 = time.time()
        while v.is_running():
            _pose(m, d, *_traj(((time.time() - t0) / period_s) % 1.0))
            v.sync()
            time.sleep(dt)


def render(n_frames=60, height=720, width=960):
    import os
    os.environ.setdefault("MUJOCO_GL", "egl")
    m, d = _model()
    try:
        r = mujoco.Renderer(m, height=height, width=width)
    except Exception:
        os.environ["MUJOCO_GL"] = "osmesa"
        r = mujoco.Renderer(m, height=height, width=width)
    imgs, poses = [], []
    for i in range(n_frames):
        ph, th = _traj(i / n_frames)
        _pose(m, d, ph, th)
        e_off = P.cell.pole_pairs * (_q(m, d, "rB") + P.r_r * th)   # B elec. vs A (lead-swapped)
        poses.append((ph, th, e_off))
        r.update_scene(d, camera="iso")
        imgs.append(r.render().copy())
    OUT.mkdir(exist_ok=True)
    from PIL import Image
    fr = [Image.fromarray(f) for f in imgs]
    fr[0].save(OUT / "sim.gif", save_all=True, append_images=fr[1:], duration=66, loop=0)
    print(f"  wrote {OUT/'sim.gif'}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    picks = [0, n_frames // 8, n_frames // 4, 3 * n_frames // 8, n_frames // 2, 5 * n_frames // 8]
    fig, ax = plt.subplots(2, 3, figsize=(12, 8))
    for a_, k in zip(ax.ravel(), picks):
        ph, th, e_off = poses[k]
        a_.imshow(imgs[k]); a_.set_axis_off()
        a_.set_title(f"hinge φ={np.degrees(ph):5.1f}°   rotor θ={np.degrees(th)%360:5.1f}°   "
                     f"B elec. offset {np.degrees(e_off):.1f}°", fontsize=9)
    fig.suptitle("2-cell Variable Chain Motor — B orbits A (stator gears) while both rotors stay "
                 "synchronized (rotor gears)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(OUT / "sim_montage.png", dpi=110)
    print(f"  wrote {OUT/'sim_montage.png'}")


if __name__ == "__main__":
    args = sys.argv[1:]
    if "check" in args:
        sys.exit(0 if check() else 1)
    elif "render" in args:
        render()
    else:
        view()
