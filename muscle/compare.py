"""
SMA cell vs HASEL cell on the same control test (option A in README.md).

Both are sized as the same 10 N / 10 mm unit and put through slowmanifold.py's test:
learn a static map from slow forced data, then track a reference with PI,
feed-forward, feed-forward + PI and lead + PI. Output is normalized 0..1 for both
(force for the SMA strand, stroke against a spring load for the HASEL unit).

Three comparisons:
  1. equal slowness rho: the paper's framing. Each cell gets a reference scaled to
     its own time constant.
  2. equal real time: the same reference, in seconds, for both. This is the design
     question.
  3. what decides whether feed-forward helps: the HASEL's time constant and the
     feedback delay, both [ASSUMED], varied.

    python muscle/compare.py        (or: make muscle-compare; takes a few minutes)
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import hasel  # noqa: E402
import slowmanifold as sm  # noqa: E402

HASEL_LOOP = 0.02                    # play half-width, fraction of v_max [ASSUMED]
PERIODS = (113.0, 34.0, 14.0, 3.4, 1.0)   # s, of a sinusoid with the same rho
PAPER = {"PI": 7.63, "FF": 3.63, "FF+PI": 2.38}   # deg RMS, Fig. 8


def cells():
    return {"SMA": sm.HystCell(), "HASEL": hasel.plant(loop=HASEL_LOOP)}


def row(r, rho):
    return " ".join(f"{100 * r['rms'][rho, c]:8.2f}" for c in sm.CONTROLLERS)


def best(r, rho):
    c = min(sm.CONTROLLERS, key=lambda c: r["rms"][rho, c])
    return c, 100 * r["rms"][rho, c]


def main():
    head = " ".join(f"{c:>8}" for c in sm.CONTROLLERS)

    print("=== 1. equal slowness: tracking RMS error, % of full output ===")
    same_rho = {}
    for name, cell in cells().items():
        r = same_rho[name] = sm.experiment(cell)
        print(f"\n  {name} (tau {cell.tau:g} s)   map error "
              + "  ".join(f"rho {rho}: {100 * r['pred'][rho]:.1f}%" for rho in sm.RHOS))
        print(f"  {'rho':>5} {'period':>8} {head}")
        for rho in sm.RHOS:
            print(f"  {rho:5.2f} {2 * np.pi * cell.tau / rho:7.1f}s {row(r, rho)}")
    h = same_rho["HASEL"]["rms"]
    print("\n  HASEL at rho 0.15 against the paper's Fig. 8, each relative to FF+PI:")
    print("    " + "   ".join(f"{c}: ours {h[0.15, c] / h[0.15, 'FF+PI']:.1f}x, paper "
                              f"{PAPER[c] / PAPER['FF+PI']:.1f}x" for c in ("PI", "FF")))

    print("\n=== 2. equal real time: the same reference for both ===")
    print(f"  {'period':>7} | {'SMA rho':>8} {'best':>8} {'error':>7} | {'HASEL rho':>9} {'best':>8} {'error':>7}")
    real = {}
    for name, cell in cells().items():
        rhos = tuple(2 * np.pi * cell.tau / p for p in PERIODS)
        real[name] = (rhos, sm.experiment(cell, rhos))
    for i, p in enumerate(PERIODS):
        out = f"  {p:6.1f}s |"
        for name in ("SMA", "HASEL"):
            rhos, r = real[name]
            c, e = best(r, rhos[i])
            out += f" {rhos[i]:{8 if name == 'SMA' else 9}.3f} {c:>8} {e:6.2f}% |"
        print(out)

    print("\n=== 3. when does feed-forward help? (rho 0.15; PI vs FF+PI, % of full output) ===")
    print(f"  {'case':38} {'PI':>7} {'FF+PI':>7}")
    cases = [
        ("HASEL tau 50 ms, delay 100 ms", hasel.plant(loop=HASEL_LOOP), sm.DELAY),
        ("HASEL tau 50 ms, delay 10 ms", hasel.plant(loop=HASEL_LOOP), 0.01),
        ("HASEL tau 200 ms, delay 100 ms", hasel.plant(loop=HASEL_LOOP, tau=0.2), sm.DELAY),
        ("HASEL tau 20 ms, delay 100 ms", hasel.plant(loop=HASEL_LOOP, tau=0.02), sm.DELAY),
        ("HASEL no hysteresis, delay 100 ms", hasel.plant(loop=0.0), sm.DELAY),
        ("SMA tau 2.7 s, delay 100 ms", sm.HystCell(), sm.DELAY),
        ("SMA no hysteresis, delay 100 ms", sm.HystCell(loop=0.0), sm.DELAY),
    ]
    for label, cell, delay in cases:
        r = sm.experiment(cell, (0.15,), delay)
        print(f"  {label:38} {100 * r['rms'][0.15, 'PI']:7.2f} {100 * r['rms'][0.15, 'FF+PI']:7.2f}")

    plot(same_rho, real)
    print(f"\nwrote {sm.OUT / 'compare.png'}")


def plot(same_rho, real):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for (name, r), col in zip(same_rho.items(), ("tab:red", "tab:blue")):
        m = r["model"]
        lo, hi = r["cell"].u_range()
        ax[0].plot((m["u"][::40] - lo) / (hi - lo), m["a"][::40], ".", ms=2, color=col, alpha=0.3)
        ax[0].plot((m["grid"][0] - lo) / (hi - lo), m["grid"][1], color=col, lw=1.8, label=f"{name}: fitted map")
    ax[0].set(xlabel="input, fraction of its working span", ylabel="output (normalized)",
              title="Static map fitted from slow forced data", ylim=(-0.1, 1.1))
    ax[0].legend(fontsize=8)

    i = PERIODS.index(3.4)
    for name, col in (("SMA", "tab:red"), ("HASEL", "tab:blue")):
        rhos, r = real[name]
        c, _ = best(r, rhos[i])
        ref, a, _ = r["trace"][rhos[i], c]
        t = np.arange(len(ref)) * r["cell"].dt
        if name == "SMA":
            ax[1].plot(t, ref, "k--", lw=1, label="reference")
        ax[1].plot(t, a, color=col, lw=1.2, label=f"{name}, {c}")
    ax[1].set(xlabel="time (s)", ylabel="output (normalized)", xlim=(t[len(t) // 3], t[-1]),
              title="The same reference for both (about 3 s per swing)")
    ax[1].legend(fontsize=8)

    w = 0.38
    for k, (name, col) in enumerate((("SMA", "tab:red"), ("HASEL", "tab:blue"))):
        rhos, r = real[name]
        ax[2].bar(np.arange(len(PERIODS)) + (k - 0.5) * w, [best(r, x)[1] for x in rhos], w,
                  color=col, label=f"{name}, best controller")
    ax[2].set(xticks=range(len(PERIODS)), xticklabels=[f"{p:g} s" for p in PERIODS],
              xlabel="reference period", ylabel="RMS error (% of full output)",
              title="Tracking error at equal real time")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    sm.OUT.mkdir(exist_ok=True)
    fig.savefig(sm.OUT / "compare.png", dpi=120)


if __name__ == "__main__":
    main()
