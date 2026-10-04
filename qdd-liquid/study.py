"""
What limits the liquid-cooled motor, with the field solved (fea.py) instead of an
assumed saturation knee, and what trading copper for iron buys at a fixed size.

For each variant of design.Design:
  torque(current)   from the nonlinear field solution
  heat(current)     from design.py's copper, drive, pump and coolant model
and so torque against heat, the current the cooling can carry, and where torque per
amp has fallen off.

    python qdd-liquid/study.py        (or: make qdd-liquid-fea; a few minutes)
"""

import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import PolyCollection  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402
import fea  # noqa: E402

J_LEVELS = np.array([3, 8, 14, 20, 26, 33, 42, 55, 75, 100, 140]) * 1e6   # A/m2 rms
HEATS = (300, 850, 1500)                                                  # W


def pick_turns(d, j=26e6):
    """Turns per coil with the least heat at a reference current (bore size against
    drive loss), among those the flow cap can cool."""
    best = None
    for n in (1, 2, 3, 4, 6, 8):
        r = dz.cooled(replace(d, turns=n), j)
        if r is not None and (best is None or dz.total_loss(r) < dz.total_loss(best[1])):
            best = (n, r)
    return replace(d, turns=best[0]) if best else d


class Variant:
    def __init__(self, label, d, steel="silicon steel", **motor_kw):
        self.label, self.d, self.steel = label, pick_turns(d), steel
        self.motor = fea.Motor(self.d, steel=steel, **motor_kw)
        self.torque = self.motor.curve([fea.amp_turns_of(self.d, j) for j in J_LEVELS])
        self.kt0 = self.torque[0] / J_LEVELS[0]                           # N·m per A/m2, unsaturated
        # heat against current, within what the capped pump can cool
        self.js = np.geomspace(2e6, 140e6, 70)
        self.cool = [dz.cooled(self.d, j) for j in self.js]
        ok = [c is not None for c in self.cool]
        self.j_cool = self.js[ok][-1]
        self.heat = np.array([dz.total_loss(c) for c in self.cool if c is not None])
        self.j_ok = self.js[ok]

    def t_at_j(self, j):
        return float(np.interp(j, J_LEVELS, self.torque))

    def kt_ratio(self, j):
        return self.t_at_j(j) / (self.kt0 * j)

    def at_heat(self, watts):
        """(torque, current density) at a heat budget, or None if cooling can't reach it."""
        if watts > self.heat[-1]:
            return None
        j = float(np.interp(watts, self.heat, self.j_ok))
        return self.t_at_j(j), j

    def j_kt_drop(self, drop=0.2):
        js = np.geomspace(J_LEVELS[0], J_LEVELS[-1], 400)
        k = np.array([self.kt_ratio(j) for j in js])
        return float(js[np.argmax(k < 1 - drop)]) if (k < 1 - drop).any() else None


def field_plot(ax, motor, amp_turns, gamma, title):
    s = motor.solve(amp_turns, gamma)
    xy, tri = motor.xy * 1e3, motor.tri
    b = np.where(motor.iron, s["b"], np.nan)
    pc = PolyCollection(xy[tri], array=b, cmap="inferno", edgecolors="none")
    pc.set_clim(0, 2.4)
    ax.add_collection(pc)
    for reg, col in ((fea.COIL, "#c87a2f"), (fea.MAGNET, "#7aa6c2")):
        ax.add_collection(PolyCollection(xy[tri[motor.region == reg]], facecolors=col, edgecolors="none"))
    r = motor.r["outer"]
    ax.set(xlim=(r * 0.42, r * 1.01), ylim=(-r * 0.02, r * 0.5), aspect="equal", title=title)
    ax.axis("off")
    return pc, s


def main():
    base = dz.Design()
    variants = [
        Variant("baseline: slots 50% of pitch, 14 mm deep", base),
        Variant("slots 40% of pitch", replace(base, slot_frac=0.4)),
        Variant("slots 30% of pitch", replace(base, slot_frac=0.3)),
        Variant("slots 10 mm deep (bigger rotor)", replace(base, slot_depth=10.0)),
        Variant("slots 40%, 10 mm deep", replace(base, slot_frac=0.4, slot_depth=10.0)),
        Variant("slots 30%, 7 mm deep", replace(base, slot_frac=0.3, slot_depth=7.0)),
        Variant("yoke 8 mm (slots 11 mm)", replace(base, stator_yoke=8.0, slot_depth=11.0)),
        Variant("magnets 5 mm", replace(base, magnet=5.0)),
        Variant("cobalt-iron laminations", base, steel="cobalt iron"),
        Variant("outrunner", replace(base, outrunner=True)),
        Variant("outrunner, slots 40%, 10 mm deep", replace(base, outrunner=True, slot_frac=0.4, slot_depth=10.0)),
    ]
    b = variants[0]

    print("=== baseline: the field solution against design.py's assumed knee ===")
    print(f"  {'A/mm²':>6} {'FEA N·m':>8} {'per amp':>8} {'knee model':>11} {'heat':>7}")
    for j in J_LEVELS:
        c = dz.cooled(b.d, j)
        knee = dz.evaluate(b.d, j)["torque"]
        heat = f"{dz.total_loss(c):6.0f}W" if c else "  (beyond the 3 L/min pump)"
        print(f"  {j/1e6:6.0f} {b.t_at_j(j):8.1f} {b.kt_ratio(j)*100:7.0f}% {knee:11.1f} {heat}")
    jd = b.j_kt_drop()
    print(f"  torque per amp is down 20% at {jd/1e6:.0f} A/mm² ({b.t_at_j(jd):.1f} N·m); "
          f"the pump's flow cap is reached at {b.j_cool/1e6:.0f} A/mm² ({b.t_at_j(b.j_cool):.1f} N·m)")

    print("\n=== same Ø120 x 25 mm stack, copper traded for iron ===")
    print(f"  {'variant':40} {'turns':>5} {'gap r':>6} {'slot':>6} " + " ".join(f"{f'{h} W':>8}" for h in HEATS)
          + f" {'cooling limit':>20} {'Kt -20% at':>11}")
    for v in variants:
        cells = []
        for h in HEATS:
            r = v.at_heat(h)
            cells.append(f"{r[0]:8.1f}" if r else f"{'-':>8}")
        jd = v.j_kt_drop()
        print(f"  {v.label:40} {v.d.turns:5d} {v.d.r_gap*1e3:5.1f}  {v.d.slot_area*1e6:5.0f}  " + " ".join(cells)
              + f" {v.t_at_j(v.j_cool):6.1f} N·m @ {v.j_cool/1e6:3.0f} A/mm²"
              + (f" {jd/1e6:7.0f} A/mm²" if jd else f" {'>140':>11}"))
    print("  (columns: torque in N·m at that total heat; gap r in mm; slot area in mm²)")

    best = max(variants[:9], key=lambda v: v.at_heat(850)[0] if v.at_heat(850) else 0)
    print(f"\n  best stator-outside variant at 850 W: {best.label}: {best.at_heat(850)[0]:.1f} N·m "
          f"against the baseline's {b.at_heat(850)[0]:.1f}")

    fig = plt.figure(figsize=(15, 8.4))
    ax = fig.add_subplot(2, 3, 1)
    ax.plot(J_LEVELS / 1e6, b.torque, "o-", lw=2, label="field solution")
    ax.plot(J_LEVELS / 1e6, b.kt0 * J_LEVELS, "--", color="gray", label="no saturation")
    ax.plot(J_LEVELS / 1e6, [dz.evaluate(b.d, j)["torque"] for j in J_LEVELS], ":", color="tab:red", label="assumed knee (design.py)")
    ax.axvline(b.j_cool / 1e6, color="tab:green", lw=1, label="3 L/min cooling limit")
    ax.set(xlabel="copper current density (A/mm² rms)", ylabel="motor torque (N·m)", ylim=(0, 45),
           title="Baseline: torque against current")
    ax.legend(fontsize=8)

    ax = fig.add_subplot(2, 3, 2)
    for v in (variants[0], variants[1], variants[4], variants[5], variants[8], variants[9]):
        ax.plot(v.heat, [v.t_at_j(j) for j in v.j_ok], lw=1.8, label=v.label)
    ax.set(xlabel="total heat (W)", ylabel="motor torque (N·m)", title="Torque against heat, to each cooling limit")
    ax.legend(fontsize=7)

    ax = fig.add_subplot(2, 3, 3)
    for v in (variants[0], variants[2], variants[5], variants[8]):
        ax.plot(J_LEVELS / 1e6, [v.kt_ratio(j) * 100 for j in J_LEVELS], "o-", ms=3, label=v.label)
    ax.set(xlabel="copper current density (A/mm² rms)", ylabel="torque per amp (% of unsaturated)",
           title="How fast each one saturates")
    ax.legend(fontsize=7)

    g, _ = b.motor.best_gamma(fea.amp_turns_of(b.d, 3e6))
    for k, j in enumerate((0.0, 26e6, 100e6)):
        ax = fig.add_subplot(2, 3, 4 + k)
        pc, s = field_plot(ax, b.motor, fea.amp_turns_of(b.d, j), g,
                           f"{j/1e6:.0f} A/mm²: iron flux density")
    fig.colorbar(pc, ax=fig.axes[3:], shrink=0.7, label="|B| in iron (T)")
    (HERE / "out").mkdir(exist_ok=True)
    fig.savefig(HERE / "out" / "study.png", dpi=110, bbox_inches="tight")
    print(f"\nwrote {HERE / 'out' / 'study.png'}")


if __name__ == "__main__":
    main()
