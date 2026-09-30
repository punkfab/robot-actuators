"""
Fit real2sim corrections on top of BAM's STS3215 model from recorded logs.

BAM's M6 was identified on Rhoban's own servos, on a pendulum. Yours differ (unit
to unit, wear, supply, the arm's own links), so fit three multipliers on the
M6 values, the ones Microduck randomizes over:

  friction_scale   velocity-independent friction budget
  gain_scale       firmware P path (error_gain_ratio): kp, duty mapping, back-EMF
  armature_scale   reflected rotor inertia

Nelder–Mead on the log-multipliers, minimizing mean joint MAE over the training
logs; the last log (or --val) is held out, as BAM does with its 75/25 split.

    python sts3215/fit.py ~/so101_logs/*.npz              # -> out/fit.json
    python sts3215/fit.py --synthetic                     # should recover ×1.4 / ×0.85
"""

import json
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from replay import SYNTH_TRUTH, load_log, mae_deg, make_synthetic, simulate  # noqa: E402

OUT = HERE / "out"
NAMES = ("friction_scale", "gain_scale", "armature_scale")


def score(logs, **kw):
    return float(np.mean([mae_deg(l, simulate(l, "m6", **kw)) for l in logs]))


def main():
    synthetic = "--synthetic" in sys.argv
    paths = make_synthetic() if synthetic else [Path(p) for p in sys.argv[1:] if not p.startswith("--")]
    if len(paths) < 2:
        print(__doc__)
        return
    logs = [load_log(p) for p in paths]
    train, val = logs[:-1], logs[-1:]
    print(f"train: {', '.join(l['name'] for l in train)}\nval:   {val[0]['name']}")

    evals = []

    def obj(x):
        kw = dict(zip(NAMES, np.exp(x)))
        s = score(train, **kw)
        evals.append(s)
        if len(evals) % 10 == 0:
            print(f"  eval {len(evals):3d}  MAE {s:.3f}°  " +
                  "  ".join(f"{n.split('_')[0]} ×{v:.3f}" for n, v in kw.items()))
        return s

    base_train, base_val = score(train), score(val)
    res = minimize(obj, np.zeros(len(NAMES)), method="Nelder-Mead",
                   options=dict(xatol=0.01, fatol=0.002, maxfev=120, initial_simplex=
                                np.vstack([np.zeros(3), 0.3 * np.eye(3)])))
    fitted = dict(zip(NAMES, map(float, np.exp(res.x))))
    fit_train, fit_val = score(train, **fitted), score(val, **fitted)

    print("\n=== fit (BAM STS3215 M6 × corrections) ===")
    for n in NAMES:
        truth = f"   truth ×{SYNTH_TRUTH[n]:.2f}" if synthetic else ""
        print(f"  {n:15s} ×{fitted[n]:.3f}{truth}")
    print(f"  train MAE {base_train:.2f}° -> {fit_train:.2f}°   "
          f"val MAE {base_val:.2f}° -> {fit_val:.2f}°   ({len(evals)} evals)")
    OUT.mkdir(exist_ok=True)
    name = "fit_synthetic.json" if synthetic else "fit.json"
    (OUT / name).write_text(json.dumps(dict(model="m6", **fitted,
        train=[l["name"] for l in train], val=[l["name"] for l in val],
        mae_train=fit_train, mae_val=fit_val), indent=2))
    print(f"wrote {OUT/name}  (use: Arm('m6', **json) minus the bookkeeping keys)")


if __name__ == "__main__":
    main()
