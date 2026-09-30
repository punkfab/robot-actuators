"""
Replay recorded SO-101 trajectories in MuJoCo and score each actuator model.

    python sts3215/replay.py ~/so101_logs/*.npz          # real logs from record.py
    python sts3215/replay.py --synthetic                  # pipeline self-test

Each log's goals are played at the sim timestep (zero-order hold, as the real bus
holds a goal until the next write); the other joints hold the recorded base pose. The
score is the mean absolute error (°) of the moving joint against the measurement,
the same metric BAM identifies on.

--synthetic makes stand-in "real" logs from BAM M6 with deliberately wrong physics
(friction ×1.4, P gain ×0.85, 4096-count encoder quantization). Replay should then
rank M6 best but not perfect, and fit.py should recover the ×1.4 / ×0.85.
"""

import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arm import ARM, TIMESTEP, Arm  # noqa: E402

OUT = HERE / "out"
SYNTH = OUT / "synthetic_logs"
SYNTH_TRUTH = dict(friction_scale=1.4, gain_scale=0.85, armature_scale=1.0)
TICK_DEG = 360.0 / 4096


def load_log(path):
    z = np.load(path, allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    d = z["data"]
    extra = z["extra"] if "extra" in z.files else None
    vin = None
    if extra is not None and extra.size and np.isfinite(extra[:, 2]).any():
        vin = float(np.nanmedian(extra[:, 2])) / 10.0
    return dict(name=Path(path).stem, t=d[:, 0], goal=np.radians(d[:, 1]),
                q=np.radians(d[:, 2:2 + len(ARM)]), joint=meta["joint"],
                base={j: np.radians(v) for j, v in meta["base_pose_deg"].items()},
                payload=meta.get("payload_g", 0.0) / 1000.0, vin=vin)


def simulate(log, kind, **kw):
    """Sim trace of all arm joints at the log's sample times."""
    if log["vin"] and kind != "stock":
        kw.setdefault("vin", log["vin"])
    a = Arm(kind, payload_kg=log["payload"], **kw)
    start = {j: log["q"][0, i] for i, j in enumerate(ARM)}
    a.reset(start)
    a.target.update(log["base"])
    t, goal, j = log["t"], log["goal"], log["joint"]
    out = np.empty_like(log["q"])
    k = 0
    for i, ts in enumerate(t):
        while k * TIMESTEP < ts:
            gi = max(0, np.searchsorted(t, k * TIMESTEP, side="right") - 1)
            a.step({j: goal[gi]})
            k += 1
        out[i] = a.q[:len(ARM)]
    return out


def mae_deg(log, sim):
    i = ARM.index(log["joint"])
    return float(np.degrees(np.mean(np.abs(sim[:, i] - log["q"][:, i]))))


def make_synthetic():
    """Stand-in real logs: M6 with wrong physics + encoder quantization."""
    sys.path.insert(0, str(HERE))
    from record import BASE_POSE_DEG, RATE_HZ, trajectories
    SYNTH.mkdir(parents=True, exist_ok=True)
    paths = []
    for traj in ("chirp", "ramp", "steps", "sub"):
        f, T = trajectories(traj)
        t = np.arange(0, T, 1 / RATE_HZ)
        goal = np.array([BASE_POSE_DEG["elbow_flex"] + f(x) for x in t])
        log = dict(t=t, goal=np.radians(goal), joint="elbow_flex",
                   base={j: np.radians(v) for j, v in BASE_POSE_DEG.items()},
                   payload=0.1, vin=7.4,
                   q=np.tile(np.radians([BASE_POSE_DEG[j] for j in ARM]), (len(t), 1)))
        q = simulate(log, "m6", **SYNTH_TRUTH)
        q = np.round(np.degrees(q) / TICK_DEG) * TICK_DEG          # encoder counts
        meta = dict(joint="elbow_flex", traj=traj, payload_g=100.0, kp=16,
                    base_pose_deg=BASE_POSE_DEG, joints=list(ARM), rate_hz=RATE_HZ,
                    synthetic_truth=SYNTH_TRUTH)
        rows = np.column_stack([t, goal, q])
        extra = np.column_stack([np.full(len(t), np.nan)] * 2 + [np.full(len(t), 74.0)])
        p = SYNTH / f"elbow_flex_{traj}_100g_synthetic.npz"
        np.savez(p, data=rows, extra=extra, meta=json.dumps(meta))
        paths.append(p)
    return paths


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    paths = make_synthetic() if "--synthetic" in sys.argv else [Path(p) for p in args]
    if not paths:
        print(__doc__)
        return
    logs = [load_log(p) for p in paths]
    kinds = ("stock", "m1", "m6")
    print(f"  {'log':34s} " + " ".join(f"{k:>7s}" for k in kinds) + "   (MAE °, moving joint)")
    sims, table = {}, []
    for log in logs:
        row = []
        for k in kinds:
            s = simulate(log, k)
            sims[(log["name"], k)] = s
            row.append(mae_deg(log, s))
        table.append(row)
        print(f"  {log['name'][:34]:34s} " + " ".join(f"{v:7.2f}" for v in row))
    mean = np.mean(table, axis=0)
    print(f"  {'mean':34s} " + " ".join(f"{v:7.2f}" for v in mean))
    plot(logs, sims, kinds)
    print(f"wrote {OUT/'replay.png'}")


def plot(logs, sims, kinds):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(len(logs), 1, figsize=(10, 2.4 * len(logs)), squeeze=False)
    for ax, log in zip(axes[:, 0], logs):
        i = ARM.index(log["joint"])
        ax.plot(log["t"], np.degrees(log["goal"]), "k:", lw=0.8, label="goal")
        ax.plot(log["t"], np.degrees(log["q"][:, i]), "k", lw=1.6, label="measured")
        for k in kinds:
            ax.plot(log["t"], np.degrees(sims[(log["name"], k)][:, i]), lw=1, label=k)
        ax.set_title(log["name"], fontsize=9)
        ax.set_ylabel(f"{log['joint']} (°)")
        ax.grid(alpha=0.3)
    axes[0, 0].legend(ncol=5, fontsize=8)
    axes[-1, 0].set_xlabel("time (s)")
    fig.tight_layout()
    fig.savefig(OUT / "replay.png", dpi=110)


if __name__ == "__main__":
    main()
