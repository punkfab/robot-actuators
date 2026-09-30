"""
How much does the STS3215 actuator model change what an SO-101 sim predicts?

stock (upstream MJCF ideal servo) vs BAM M1 (firmware P law + Coulomb/viscous) vs
BAM M6 (+ Stribeck, load-dependent, directional, quadratic friction). All BAM runs
use lerobot's firmware P = 16 at 7.4 V.

  A. static sag      hold a pose with 0 / 100 / 200 g at the gripper
  B. hysteresis      arrive at the same pose from above and from below
  C. step response   0.5 rad elbow step: rise, overshoot, final error
  D. uncertainty     M6 with friction ×[0.7, 1.3], supply [6.8, 8.0] V (Microduck-style DR)
  E. creep           M6 held 10 s with MuJoCo's default soft friction vs stiffened

  python sts3215/compare.py        # report + out/compare.png
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arm import ARM, POSES, TIMESTEP, Arm  # noqa: E402

OUT = HERE / "out"
MODELS = ("stock", "m1", "m6")


def tip_err_mm(a, pose):
    return float(np.linalg.norm(a.tip() - a.tip_at(pose)) * 1000)


def sag(kind, pose_name, payload, seconds=3.0, **kw):
    a = Arm(kind, payload_kg=payload, **kw)
    pose = POSES[pose_name]
    a.reset(pose)
    a.run(seconds)
    dq = np.degrees(a.q[:5] - np.array([pose[j] for j in ARM]))
    return tip_err_mm(a, pose), dq


def approach(kind, pose_name, sign, payload=0.1, offset=0.25):
    """Start offset (rad) above/below on shoulder+elbow, command the pose, settle."""
    pose = POSES[pose_name]
    a = Arm(kind, payload_kg=payload)
    start = dict(pose)
    start["shoulder_lift"] += sign * offset
    start["elbow_flex"] -= sign * offset
    a.reset(start)
    a.run(0.5)
    a.run(2.5, target=pose)
    return a.tip() - a.tip_at(pose)


def step(kind, joint="elbow_flex", size=0.5, seconds=1.0):
    a = Arm(kind, payload_kg=0.1)
    pose = POSES["mid"]
    a.reset(pose)
    a.run(0.5)
    i = ARM.index(joint)
    q0 = a.q[i]
    tgt = pose[joint] + size
    t, y = [], []
    n = int(seconds / TIMESTEP)
    for k in range(n):
        a.step({joint: tgt} if k == 0 else None)
        t.append(k * TIMESTEP)
        y.append(a.q[i])
    t, y = np.array(t), np.array(y)
    frac = (y - q0) / (tgt - q0)
    q0 = pose[joint]                       # plot from the commanded start, not the sagged one
    rise = t[np.argmax(frac >= 0.9)] - t[np.argmax(frac >= 0.1)] if frac.max() >= 0.9 else np.nan
    return t, y - q0, dict(rise_ms=rise * 1000, overshoot=max(0.0, frac.max() - 1) * 100,
                          final_err_deg=np.degrees(tgt - y[-1]))


def main():
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(0)

    print("=== A. static sag: gripper-frame error after 3 s holding (mm) ===")
    print(f"  {'pose':6s} {'payload':>7s}  " + "  ".join(f"{m:>7s}" for m in MODELS)
          + "   M6 joint error (lift/elbow/wrist °)")
    sag_rows = {}
    for pose in ("reach", "mid"):
        for pay in (0.0, 0.1, 0.2):
            errs = []
            for kind in MODELS:
                e, dq = sag(kind, pose, pay)
                errs.append(e)
            sag_rows[(pose, pay)] = errs
            print(f"  {pose:6s} {pay*1000:5.0f} g  " + "  ".join(f"{e:7.1f}" for e in errs)
                  + f"   {dq[1]:+.2f} / {dq[2]:+.2f} / {dq[3]:+.2f}")
    print()

    print("=== B. hysteresis: arrive at 'reach' (100 g) from above vs below ===")
    for kind in MODELS:
        up, down = approach(kind, "reach", +1), approach(kind, "reach", -1)
        print(f"  {kind:6s} from above {np.linalg.norm(up)*1000:5.1f} mm  from below "
              f"{np.linalg.norm(down)*1000:5.1f} mm  -> band {np.linalg.norm(up-down)*1000:5.1f} mm")
    print("  (a stock sim returns to the same point from either side; friction makes the"
          " rest point depend on the approach — that band is the repeatability floor)\n")

    print("=== C. step response: elbow +0.5 rad at 'mid', 100 g ===")
    traces = {}
    for kind in MODELS:
        t, y, s = step(kind)
        traces[kind] = (t, y)
        print(f"  {kind:6s} rise(10–90%) {s['rise_ms']:6.0f} ms  overshoot {s['overshoot']:4.1f} %  "
              f"final error {s['final_err_deg']:+.2f}°")
    print("  (BAM's firmware rate-limits the internal target to ~5.2 rad/s, so big steps ramp)\n")

    print("=== D. M6 uncertainty band (friction ×0.7–1.3, supply 6.8–8.0 V), reach + 100 g ===")
    band = []
    for _ in range(24):
        fs, v = rng.uniform(0.7, 1.3), rng.uniform(6.8, 8.0)
        band.append(sag("m6", "reach", 0.1, friction_scale=fs, vin=v)[0])
    band = np.array(band)
    print(f"  tip error {band.min():.1f}–{band.max():.1f} mm (median {np.median(band):.1f}); "
          f"stock says {sag_rows[('reach', 0.1)][0]:.1f} mm\n")

    print("=== E. creep: M6, reach + 200 g, 10 s hold ===")
    for stiff in (False, True):
        e3, _ = sag("m6", "reach", 0.2, seconds=3.0, stiff_friction=stiff)
        e10, _ = sag("m6", "reach", 0.2, seconds=10.0, stiff_friction=stiff)
        label = "stiffened solimp" if stiff else "MuJoCo default   "
        print(f"  {label}  3 s {e3:5.1f} mm   10 s {e10:5.1f} mm   drift {e10-e3:+.1f} mm")
    print()

    plot(sag_rows, traces, band)
    print(f"wrote {OUT/'compare.png'}")


def plot(sag_rows, traces, band):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 4))
    labels = [f"{p}\n{int(w*1000)} g" for p, w in sag_rows]
    x = np.arange(len(labels))
    for k, kind in enumerate(MODELS):
        ax1.bar(x + (k - 1) * 0.27, [v[k] for v in sag_rows.values()], 0.27,
                label={"stock": "stock MJCF (ideal servo)", "m1": "BAM M1",
                       "m6": "BAM M6"}[kind])
    i = labels.index("reach\n100 g")
    ax1.errorbar(i + 0.27, np.median(band), yerr=[[np.median(band) - band.min()],
                 [band.max() - np.median(band)]], fmt="none", ecolor="k", capsize=3,
                 label="M6 friction/voltage band")
    ax1.set_xticks(x, labels)
    ax1.set_ylabel("gripper sag after 3 s hold (mm)")
    ax1.set_title("Holding a pose: what each actuator model predicts")
    ax1.legend(fontsize=8)
    ax1.grid(axis="y", alpha=0.3)
    for kind, (t, y) in traces.items():
        ax2.plot(t * 1000, np.degrees(y), label=kind)
    ax2.axhline(np.degrees(0.5), color="k", ls="--", lw=0.8)
    ax2.set_xlabel("time (ms)")
    ax2.set_ylabel("elbow displacement (°)")
    ax2.set_title("0.5 rad elbow step (mid pose, 100 g)")
    ax2.legend()
    ax2.grid(alpha=0.3)
    fig.suptitle("SO-101 / STS3215: stock MuJoCo servo vs BAM actuator models (firmware P = 16)")
    fig.tight_layout()
    fig.savefig(OUT / "compare.png", dpi=130)


if __name__ == "__main__":
    main()
