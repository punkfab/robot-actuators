"""
Designing a joint set with RobStride actuators, end to end, on a small arm:

  1. size    inverse dynamics along the task -> lightest RobStride per joint that
             passes peak, RMS (thermal) and back-EMF envelope, with the actuators'
             own masses fed back into the arm (sizing.py)
  2. tune    MIT gains from the design: inverse-dynamics feed-forward carries the
             motion, kp/kd set a closed-loop bandwidth on each joint's inertia
             (kp = J·ω², kd = 2ζ·J·ω, the rule Berkeley Humanoid Lite uses), clipped
             to the model's MIT range
  3. verify  run the task through the honest driver model (servo.py) and through an
             ideal MuJoCo servo; report tip error, saturation and thermal load, and
             run the datasheet guess until its driver trips

The arm: shoulder pitch, elbow, wrist pitch in one vertical plane (a yaw base would
add a fourth joint the same way), 0.35 + 0.30 + 0.10 m aluminium-tube links, a
payload at the tool. Task: pick at low reach (0.4 s dwell), lift over, place at full
reach (0.4 s dwell), back: 3.2 s cycles.

    python robstride/demo.py            # prints the tables, writes out/demo.png
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import mujoco  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from servo import RobStrideServos, motorize  # noqa: E402
from sizing import add_actuator_masses, check_picks, size_joints  # noqa: E402
from specs import SPECS  # noqa: E402

OUT = HERE / "out"
JOINTS = ("shoulder", "elbow", "wrist")
LINKS = (0.35, 0.30, 0.10)          # m
LINK_KG_PER_M = 0.9                 # 30×2 mm Al tube + clamps, roughly
PAYLOAD = 2.0                       # kg at the tool
DT = 0.001
MOVE = 0.6                          # s per segment
DWELL = 0.4                         # s held at pick and at place
CYCLE = 4 * MOVE + 2 * DWELL

# Waypoints (shoulder, elbow, wrist), rad; shoulder 0 = arm horizontal, +up.
PICK = np.array([-0.9, 1.6, 0.9])
OVER = np.array([0.6, 0.9, -0.3])
PLACE = np.array([0.05, 0.05, -0.1])


def make_spec(payload=PAYLOAD):
    s = mujoco.MjSpec()
    s.option.timestep = DT
    s.option.integrator = mujoco.mjtIntegrator.mjINT_IMPLICITFAST
    body = s.worldbody.add_body(name="base", pos=[0, 0, 0.5])
    body.add_geom(type=mujoco.mjtGeom.mjGEOM_BOX, size=[0.05, 0.05, 0.05], mass=2.0,
                  contype=0, conaffinity=0)
    parent, x = body, 0.0
    for name, L in zip(JOINTS, LINKS):
        b = parent.add_body(name=name, pos=[x, 0, 0])
        b.add_joint(name=name, type=mujoco.mjtJoint.mjJNT_HINGE, axis=[0, -1, 0])
        b.add_geom(type=mujoco.mjtGeom.mjGEOM_CAPSULE, fromto=[0, 0, 0, L, 0, 0],
                   size=[0.015, 0, 0], mass=LINK_KG_PER_M * L, contype=0, conaffinity=0)
        parent, x = b, L
    tool = parent.add_body(name="tool", pos=[x, 0, 0])
    tool.add_site(name="tip")
    if payload:
        tool.add_geom(type=mujoco.mjtGeom.mjGEOM_SPHERE, size=[0.03, 0, 0], mass=payload,
                      contype=0, conaffinity=0, rgba=[0.9, 0.5, 0.1, 1])
    return s


def min_jerk(a, b, n):
    s = np.linspace(0, 1, n)[:, None]
    return a + (b - a) * (10 * s**3 - 15 * s**4 + 6 * s**5)


def task(cycles=1):
    n, h = int(MOVE / DT), int(DWELL / DT)
    one = np.vstack([np.repeat(PICK[None], h, 0), min_jerk(PICK, OVER, n),
                     min_jerk(OVER, PLACE, n), np.repeat(PLACE[None], h, 0),
                     min_jerk(PLACE, OVER, n), min_jerk(OVER, PICK, n)])
    return np.vstack([one] * cycles)


def build(picks):
    spec = make_spec()
    add_actuator_masses(spec, picks)
    motorize(spec, JOINTS)
    return spec.compile()


def tune_gains(picks, bw_hz=6.0, zeta=0.9):
    """kp = J·ω², kd = 2ζ·J·ω on the joint's inertia at the stretched pose (link +
    distal masses + reflected armature), clipped to the MIT range."""
    m = build(picks)
    d = mujoco.MjData(m)
    mujoco.mj_forward(m, d)
    M = np.zeros((m.nv, m.nv))
    mujoco.mj_fullM(m, M, d.qM)
    w = 2 * np.pi * bw_hz
    gains = {}
    for jn in JOINTS:
        s = SPECS[picks[jn]]
        dof = m.jnt_dofadr[m.joint(jn).id]
        J = M[dof, dof] + s.armature
        gains[jn] = (min(J * w * w, s.kp_max), min(2 * zeta * J * w, s.kd_max))
    return gains


def verify(picks, gains, ideal, ff=True, seconds=3 * CYCLE, **servo_kw):
    """Drive the task; ff = inverse-dynamics feed-forward (gravity + inertia +
    armature, no friction: the controller doesn't know that)."""
    m = build(picks)
    d = mujoco.MjData(m)
    q_ref = task(int(np.ceil(seconds / CYCLE)))[: int(seconds / DT)]
    dq_ref = np.gradient(q_ref, DT, axis=0)
    ddq_ref = np.gradient(dq_ref, DT, axis=0)
    d.qpos[:] = q_ref[0]
    mujoco.mj_forward(m, d)
    rs = RobStrideServos(m, d, picks, ideal=ideal, **servo_kw)  # sets armature
    for jn in JOINTS:
        rs.set(jn, kp=gains[jn][0], kd=gains[jn][1])
    rs.reset()
    # The host's model: rigid body + the datasheet armature, no friction or damping.
    # (mj_inverse on `m` itself would include the friction the updater just set, a
    # feed-forward that knows the true friction every step.)
    mc = build(picks)
    mc.dof_armature[:] = 0.0 if ideal else [SPECS[picks[j]].armature for j in JOINTS]
    dg = mujoco.MjData(mc)
    tip = m.site("tip").id
    log = dict(t=[], tau=[], theta=[], tip_err=[])
    for k in range(len(q_ref)):
        dg.qpos[:], dg.qvel[:], dg.qacc[:] = q_ref[k], dq_ref[k], ddq_ref[k]
        mujoco.mj_inverse(mc, dg)                   # also places the reference tip
        for i, jn in enumerate(JOINTS):
            rs.set(jn, p=q_ref[k, i], v=dq_ref[k, i],
                   tau_ff=dg.qfrc_inverse[i] if ff else 0.0)
        rs.step()
        log["t"].append(d.time)
        log["tau"].append([rs.j[j].tau for j in JOINTS])
        log["theta"].append([rs.j[j].thermal.theta for j in JOINTS])
        log["tip_err"].append(np.linalg.norm(d.site_xpos[tip] - dg.site_xpos[tip]))
    out = {k: np.array(v) for k, v in log.items()}
    out["fault"] = {j: rs.j[j].fault for j in JOINTS}
    return out


def summary(label, log, picks):
    last = log["t"] > log["t"][-1] - CYCLE
    sat = max(np.mean(np.abs(log["tau"][:, i]) > 0.98 * SPECS[picks[j]].tau_peak)
              for i, j in enumerate(JOINTS))
    faults = ", ".join(f"{j} tripped at {t:.0f} s" for j, t in log["fault"].items() if t)
    print(f"  {label:26s} tip error mean {1e3*np.mean(log['tip_err'][last]):5.1f} mm "
          f"max {1e3*np.max(log['tip_err'][last]):5.1f} mm   saturated {sat:3.0%}   "
          f"{faults or 'no trips'}")


def main():
    OUT.mkdir(exist_ok=True)
    kinds = sorted(SPECS, key=lambda k: SPECS[k].mass)
    q = task(1)
    print(f"=== 1. size: 3-DOF arm, {PAYLOAD:.0f} kg payload, {CYCLE:.1f} s pick-place cycle ===")
    picks, checks, every = size_joints(make_spec, JOINTS, q, DT, margin=1.2)
    for jn in JOINTS:
        print(f"  {jn}:")
        for c in every[jn]:
            print(f"    {c.row()}")
            if c.kind == picks[jn]:
                break
    total = sum(SPECS[k].mass for k in picks.values())
    print(f"  -> {'  '.join(f'{j}={k}' for j, k in picks.items())}   "
          f"actuators on the arm {total:.2f} kg")

    m0 = make_spec().compile()
    d0 = mujoco.MjData(m0)
    mujoco.mj_forward(m0, d0)                      # stretched horizontal
    naive = {}
    for jn in JOINTS:
        g = abs(d0.qfrc_bias[m0.jnt_dofadr[m0.joint(jn).id]])
        naive[jn] = next(k for k in kinds if SPECS[k].tau_peak >= 1.5 * g)
    print("\n  datasheet guess (peak ≥ 1.5 × static gravity at full reach, actuator masses "
          "ignored):\n    " + "  ".join(f"{j}={k}" for j, k in naive.items()))
    for jn, c in check_picks(make_spec, JOINTS, naive, q, DT).items():
        print(f"    {jn:8s} {c.row()}")

    print("\n=== 2. tune: ID feed-forward + kp/kd for 6 Hz, ζ 0.9 ===")
    gains = tune_gains(picks)
    for jn in JOINTS:
        print(f"  {jn:8s} {picks[jn]:5s} kp {gains[jn][0]:7.1f} N·m/rad   "
              f"kd {gains[jn][1]:5.2f} N·m·s/rad")

    print("\n=== 3. verify (last cycle of 3) ===")
    runs = {}
    for label, ideal, ff in (("MuJoCo ideal servo, PD", True, False),
                             ("RobStride model, PD", False, False),
                             ("MuJoCo ideal servo, PD+FF", True, True),
                             ("RobStride model, PD+FF", False, True)):
        runs[label] = verify(picks, gains, ideal, ff)
        summary(label, runs[label], picks)
    # What isn't known: gearbox friction (η is assumed), latency, host rate, bus sag.
    rng = np.random.default_rng(0)
    band = []
    for _ in range(12):
        kw = dict(friction_scale=np.exp(rng.uniform(np.log(0.5), np.log(2.0))),
                  latency=rng.uniform(0.3e-3, 2e-3), host_hz=rng.choice([250.0, 500.0, 1000.0]),
                  vbus=rng.uniform(40.0, 50.0), armature_scale=rng.uniform(1.0, 3.0))
        log = verify(picks, gains, False, True, **kw)
        last = log["t"] > log["t"][-1] - CYCLE
        band.append((1e3 * np.mean(log["tip_err"][last]), 1e3 * np.max(log["tip_err"][last]), kw))
    mx = np.array([b[1] for b in band])
    worst = band[int(np.argmax(mx))][2]
    print(f"  {'RobStride PD+FF, DR band':26s} tip error max {mx.min():.1f}–{mx.max():.1f} mm over 12 draws "
          f"(friction ×0.5–2, latency 0.3–2 ms, host 250–1000 Hz, 40–50 V, armature ×1–3)\n"
          f"  {'':26s} worst draw: friction ×{worst['friction_scale']:.2f}, "
          f"{worst['latency']*1e3:.1f} ms, {worst['host_hz']:.0f} Hz")

    gn = tune_gains(naive)
    long = verify(naive, gn, False, True, seconds=90.0)
    summary("datasheet guess, 90 s", long, naive)

    fig, ax = plt.subplots(3, 1, figsize=(9, 9))
    for label in ("MuJoCo ideal servo, PD+FF", "RobStride model, PD", "RobStride model, PD+FF"):
        ax[0].plot(runs[label]["t"], 1e3 * runs[label]["tip_err"], label=label)
    ax[0].set_ylabel("tool tip error (mm)")
    ax[0].set_xlabel("time (s)")
    ax[0].legend(fontsize=8)
    log = runs["RobStride model, PD+FF"]
    for i, jn in enumerate(JOINTS):
        s = SPECS[picks[jn]]
        ax[1].plot(log["t"], log["tau"][:, i], color=f"C{i}", label=f"{jn} ({s.name})")
        ax[1].axhline(s.tau_rated, color=f"C{i}", ls=":", lw=0.8)
    ax[1].set_ylabel("driver torque (N·m); dotted = rated")
    ax[1].set_xlabel("time (s)")
    ax[1].legend(fontsize=8)
    for i, jn in enumerate(JOINTS):
        ax[2].plot(long["t"], long["theta"][:, i], color=f"C{i}",
                   label=f"{jn} ({naive[jn]}, datasheet guess)")
    ax[2].axhline(1.0, color="k", ls="--", lw=0.8)
    ax[2].set_ylabel("thermal load θ (1 = trip)")
    ax[2].set_xlabel("time (s)")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "demo.png", dpi=120)
    print(f"\nwrote {OUT/'demo.png'}")


if __name__ == "__main__":
    main()
