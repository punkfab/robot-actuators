"""
Identify a RobStride's gearbox on a pendulum bench: record -> replay -> fit.

The numbers RobStride doesn't publish are the ones the demo shows matter most:
gearbox efficiency (load-dependent friction), drag, viscous loss, and the
reflected inertia (T-K-233's sysid found 1.4–15× the booklet's). BAM's recipe
identifies exactly these on a pendulum, so do the same: the actuator's output
axis horizontal, an arm with a mass at R, at two masses (25% and 60% of rated
torque): load-dependent friction needs load to show, and a second load to tell it
apart from constant drag.

Trajectories, all as MIT frames at the host rate (the same generator drives the
hardware in record.py and the sim here, so the fit replays exactly what was sent):

  coast   hold 90° (kp), then kp = kd = τ = 0: a free swing that decays by friction
          alone. Cleanest friction data; drive and backdrive both happen every swing.
  sweep   kp = kd = 0, τ = A·sin(chirp 0.3→1.5 Hz): the motor pumps the pendulum,
          alternately driving and being backdriven, across speeds
  raise   kp moderate, p* a slow triangle through ±80°: gravity load rises and falls
          at low speed; that's where the load-dependent friction shows (hold window)
  steps   kp/kd position steps of ±20°/±40°: armature and the loop together

Fitted: Kc (drag), Km, Ke (so η_drive and η_back separately), Kv, armature. First
an equation-error regression (everything enters the pendulum's dynamics linearly),
then Nelder–Mead on the replayed trajectories to polish it (multipliers on the
prior; Ke additively, its prior is 0), with a held-out log, as sts3215/fit.py does.
Trajectory-only fitting from the prior stalls: long swings make the angle error a
bumpy objective.

    python robstride/bench.py --synthetic RS03     # self-test: recover a known truth
    python robstride/bench.py RS03 logs/*.npz      # fit real logs from record.py
"""

import json
import sys
from pathlib import Path

import mujoco
import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from servo import RobStrideServos, gearbox_prior, motorize  # noqa: E402
from specs import SPECS  # noqa: E402

OUT = HERE / "out"
R = 0.25                  # m, pivot to mass centre
HOST_HZ = 500.0
DT = 0.001
LOAD_FRACTION = 0.5       # m·g·R as a fraction of the rated torque
SOFT_LIMIT = np.radians(120)
TRAJS = ("coast", "sweep", "raise", "steps")
# Two loads: at one, constant drag and load-proportional friction are nearly
# interchangeable (BAM records several masses for the same reason). Steps held out.
PLAN = (("coast", 0.25), ("coast", 0.6), ("sweep", 0.25), ("sweep", 0.6),
        ("raise", 0.6), ("steps", 0.4))
DURATION = {"coast": 6.0, "sweep": 8.0, "raise": 8.0, "steps": 6.0}


def bench_mass(kind, load=LOAD_FRACTION):
    return load * SPECS[kind].tau_rated / (9.81 * R)


def start_angle(traj):
    return np.pi / 2 if traj == "coast" else 0.0


def command(traj, t, q, dq, kind, load=LOAD_FRACTION):
    """(p, v, kp, kd, τ) the host sends at time t (q is 0 hanging down, + lifts).
    Closed-loop only for the safety clamp; everything else is a function of t."""
    s = SPECS[kind]
    mgr = bench_mass(kind, load) * 9.81 * R
    J = bench_mass(kind, load) * R * R
    kp_hold = min(4 * mgr, s.kp_max)                        # sags ~15° at 90°
    kd_hold = min(2 * 0.7 * np.sqrt(kp_hold * J), s.kd_max)
    if abs(q) > SOFT_LIMIT:                                 # past the soft limit: catch it
        return (np.sign(q) * SOFT_LIMIT * 0.8, 0.0, kp_hold, kd_hold, 0.0)
    if traj == "coast":
        if t < 1.0:
            return (np.pi / 2, 0.0, kp_hold, kd_hold, mgr)
        return (0.0, 0.0, 0.0, 0.0, 0.0)
    if traj == "sweep":
        f0, f1, T = 0.3, 1.5, DURATION["sweep"]
        ph = 2 * np.pi * (f0 * t + (f1 - f0) * t * t / (2 * T))
        return (0.0, 0.0, 0.0, 0.0, 0.2 * mgr * np.sin(ph) * min(t, 1.0))
    if traj == "raise":
        a, T = np.radians(80), DURATION["raise"]
        u = (t % (T / 2)) / (T / 2)
        tri = 1 - abs(2 * u - 1)                            # 0 -> 1 -> 0 every T/2
        sign = 1 if t < T / 2 else -1
        return (sign * a * tri, 0.0, kp_hold, kd_hold, 0.0)
    if traj == "steps":
        seq = np.radians([0, 20, -20, 40, -40, 0])
        p = seq[min(int(t / 1.0), len(seq) - 1)]
        return (p, 0.0, kp_hold / 2, kd_hold, 0.0)
    raise ValueError(traj)


def make_model(kind, load=LOAD_FRACTION):
    s = mujoco.MjSpec()
    s.option.timestep = DT
    s.option.integrator = mujoco.mjtIntegrator.mjINT_IMPLICITFAST
    arm = s.worldbody.add_body(name="arm")
    arm.add_joint(name="out", type=mujoco.mjtJoint.mjJNT_HINGE, axis=[0, 1, 0])
    arm.add_geom(type=mujoco.mjtGeom.mjGEOM_SPHERE, size=[0.03, 0, 0], pos=[0, 0, -R],
                 mass=bench_mass(kind, load), contype=0, conaffinity=0)
    arm.add_geom(type=mujoco.mjtGeom.mjGEOM_CAPSULE, fromto=[0, 0, 0, 0, 0, -R],
                 size=[0.008, 0, 0], mass=1e-3, contype=0, conaffinity=0)
    motorize(s, {"out": kind})
    return s.compile()


def run(kind, traj, friction=None, armature=None, commands=None, latency=0.5e-3,
        noise=False, rng=None, load=LOAD_FRACTION):
    """Simulate one trajectory. commands=None: generate them closed-loop (synthetic
    recording); otherwise replay the logged (N, 5) commands open-loop."""
    m = make_model(kind, load)
    d = mujoco.MjData(m)
    d.qpos[0] = start_angle(traj)
    mujoco.mj_forward(m, d)
    rs = RobStrideServos(m, d, {"out": kind}, host_hz=HOST_HZ, latency=latency,
                         friction=friction, armature=armature)
    n = int(DURATION[traj] * HOST_HZ) if commands is None else len(commands)
    sub = int(round(1 / (HOST_HZ * DT)))
    cmds, qs, dqs, taus = [], [], [], []
    for k in range(n):
        t = k / HOST_HZ
        q, dq = float(d.qpos[0]), float(d.qvel[0])
        if noise:                                            # 14-bit output encoder
            q = np.round(q / (2 * np.pi / 16384)) * (2 * np.pi / 16384)
            dq = dq + rng.normal(0, 0.01)
        c = command(traj, t, q, dq, kind, load) if commands is None else commands[k]
        rs.set("out", p=c[0], v=c[1], kp=c[2], kd=c[3], tau_ff=c[4])
        if k == 0:
            rs.reset()
        rs.step(sub)
        cmds.append(c)
        qs.append(q)
        dqs.append(dq)
        taus.append(rs.j["out"].tau)
    return dict(kind=kind, traj=traj, load=load, t=np.arange(n) / HOST_HZ, cmd=np.array(cmds),
                q=np.array(qs), dq=np.array(dqs), tau=np.array(taus))


def mae_deg(log, sim):
    return float(np.degrees(np.mean(np.abs(log["q"] - sim["q"]))))


# ---- fit ---------------------------------------------------------------------------
NAMES = ("drag", "km", "ke", "kv", "armature")


def params(kind, x):
    """Friction + armature from log-multipliers x on the spec prior."""
    s = SPECS[kind]
    f0 = gearbox_prior(s)
    e = np.exp(x)
    # Ke's prior is 0 (Katz), so it's fitted additively, in units of the prior's Km.
    f = f0.scaled(kc=f0.kc * e[0], km=f0.km * e[1], ke=max(0.0, f0.km * x[2]),
                  kv=max(f0.kv, 1e-3) * e[3])
    f.name = f"{kind} fitted"
    return f, s.armature * e[4]


def score(logs, x):
    f, a = params(logs[0]["kind"], x)
    return float(np.mean([mae_deg(l, run(l["kind"], l["traj"], f, a, commands=l["cmd"],
                                         load=float(l["load"]))) for l in logs]))


def regress(logs, v_min=0.3):
    """Equation-error least squares, linear in (armature, Kc, Km, Ke, Kv):
        τm + τe − m·R²·q̈ = J_a·q̈ + sgn(q̇)·(Kc + Km|τm| + Ke|τe|) + Kv·q̇
    τm is the driver's torque (on hardware: its reported Kt·iq), τe = −m·g·R·sin q.
    Uses moving samples where τm and τe oppose (BAM's budget is |Km·τm − Ke·τe|,
    which is Km|τm| + Ke|τe| only then). q̇, q̈ from a Savitzky–Golay fit of q."""
    from scipy.signal import savgol_filter
    A, b = [], []
    for l in logs:
        dt = 1 / HOST_HZ
        m = bench_mass(l["kind"], float(l["load"]))
        q = savgol_filter(l["q"], 31, 3)
        dq = savgol_filter(l["q"], 31, 3, deriv=1, delta=dt)
        ddq = savgol_filter(l["q"], 31, 3, deriv=2, delta=dt)
        tm, te = l["tau"], -m * 9.81 * R * np.sin(q)
        use = (np.abs(dq) > v_min) & (tm * te < 0)
        sg = np.sign(dq)
        A.append(np.c_[ddq, sg, sg * np.abs(tm), sg * np.abs(te), dq][use])
        b.append((tm + te - m * R * R * ddq)[use])
    A, b = np.vstack(A), np.concatenate(b)
    sol, *_ = np.linalg.lstsq(A, b, rcond=None)
    return np.maximum(sol, [1e-5, 0, 0, 0, 0]), len(b)


def x_from(kind, arm, kc, km, ke, kv):
    """Inverse of params(): physical values -> log-multipliers on the prior."""
    f0 = gearbox_prior(SPECS[kind])
    return np.array([np.log(max(kc, 1e-4) / f0.kc), np.log(max(km, 1e-4) / f0.km),
                     ke / f0.km, np.log(max(kv, 1e-5) / max(f0.kv, 1e-3)),
                     np.log(arm / SPECS[kind].armature)])


def fit(logs, maxfev=120, verbose=True):
    """Regression for the start, then Nelder–Mead on the trajectories to polish."""
    train, val = logs[:-1], logs[-1:]
    evals = []

    def obj(x):
        v = score(train, x)
        evals.append(v)
        if verbose and len(evals) % 20 == 0:
            print(f"  eval {len(evals):3d}  MAE {v:.3f}°")
        return v

    sol, n = regress(train)
    xr = x_from(train[0]["kind"], *sol)
    if verbose:
        print(f"  regression ({n} samples): {describe(train[0]['kind'], xr)}")
    x0 = np.zeros(len(NAMES))
    res = minimize(obj, xr, method="Nelder-Mead",
                   options=dict(xatol=0.005, fatol=0.001, maxfev=maxfev,
                                initial_simplex=np.vstack([xr, xr + 0.1 * np.eye(len(NAMES))])))
    return res.x, dict(train0=score(train, x0), val0=score(val, x0),
                       train_reg=score(train, xr),
                       train=score(train, res.x), val=score(val, res.x), evals=len(evals))


def save_log(path, log):
    np.savez(path, **{k: v for k, v in log.items()})


def load_log(path):
    z = np.load(path, allow_pickle=True)
    return {k: (z[k].item() if z[k].shape == () else z[k]) for k in z.files}


def describe(kind, x):
    f, a = params(kind, x)
    return (f"drag {f.kc:.3f} N·m  η_drive {f.eta_drive:.1%}  η_back {f.eta_back:.1%}  "
            f"kv {f.kv:.4f}  armature {a:.4f} kg·m² ({a / SPECS[kind].armature:.2f}× booklet)")


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    kind = args[0] if args else "RS03"
    OUT.mkdir(exist_ok=True)
    if "--synthetic" in sys.argv:
        # A plausible "real" unit: lossier than the prior, asymmetric, heavier rotor.
        s = SPECS[kind]
        f0 = gearbox_prior(s)
        truth_x = np.array([np.log(1.6), np.log(1.5), 0.5, np.log(2.0), np.log(2.2)])
        # drag ×1.6, Km ×1.5, Ke = 0.5·Km_prior, Kv ×2, armature ×2.2
        truth_f, truth_a = params(kind, truth_x)
        print(f"{kind} bench: {R} m arm, loads "
              + ", ".join(f"{bench_mass(kind, ld):.2f} kg ({ld:.0%} of rated)" for ld in (0.25, 0.4, 0.6)))
        print(f"  prior : {describe(kind, np.zeros(5))}")
        print(f"  truth : {describe(kind, truth_x)}")
        rng = np.random.default_rng(0)
        logs = [run(kind, tr, truth_f, truth_a, noise=True, rng=rng, load=ld,
                    latency=rng.uniform(0.3e-3, 1e-3)) for tr, ld in PLAN]
        for l in logs:
            save_log(OUT / f"bench_{kind}_{l['traj']}_{l['load']:.2f}_synthetic.npz", l)
    else:
        logs = [load_log(p) for p in args[1:]]
        if len(logs) < 2:
            print(__doc__)
            return
    order = sorted(logs, key=lambda l: l["traj"] == "steps")   # hold out 'steps' if present
    x, r = fit(order)
    print(f"  fitted: {describe(kind, x)}")
    print(f"  train MAE prior {r['train0']:.2f}° -> regression {r['train_reg']:.2f}° -> {r['train']:.2f}°   "
          f"held-out ({order[-1]['traj']}) {r['val0']:.2f}° -> {r['val']:.2f}°   ({r['evals']} evals)")
    f, a = params(kind, x)
    tag = "_synthetic" if "--synthetic" in sys.argv else ""
    (OUT / f"bench_{kind}{tag}.json").write_text(json.dumps(dict(
        kind=kind, kc=f.kc, km=f.km, ke=f.ke, kv=f.kv, armature=a,
        eta_drive=f.eta_drive, eta_back=f.eta_back, mae=r), indent=2))
    plot(kind, order, x, tag)


def plot(kind, logs, x, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    f, a = params(kind, x)
    fig, ax = plt.subplots(len(logs), 1, figsize=(9, 2.3 * len(logs)))
    for axi, l in zip(np.atleast_1d(ax), logs):
        prior = run(kind, l["traj"], commands=l["cmd"], load=float(l["load"]))
        fitted = run(kind, l["traj"], f, a, commands=l["cmd"], load=float(l["load"]))
        axi.plot(l["t"], np.degrees(l["q"]), "k", lw=2, alpha=0.35, label="log")
        axi.plot(l["t"], np.degrees(prior["q"]), lw=1, label=f"prior ({mae_deg(l, prior):.1f}°)")
        axi.plot(l["t"], np.degrees(fitted["q"]), lw=1, label=f"fitted ({mae_deg(l, fitted):.2f}°)")
        axi.set_ylabel(f"{l['traj']} {float(l['load']):.0%} (°)")
        axi.legend(fontsize=7, loc="upper right")
    np.atleast_1d(ax)[-1].set_xlabel("time (s)")
    fig.tight_layout()
    fig.savefig(OUT / f"bench_{kind}{tag}.png", dpi=110)
    print(f"  wrote {OUT / f'bench_{kind}{tag}.png'}")


if __name__ == "__main__":
    main()
