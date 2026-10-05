"""
What if the conductors are solid and the only cooling is the jacket around the stator?

Same iron and slots as design.SHALLOW, the same single column of four turns per coil
side, but solid rectangular conductors (the bore filled in: 27% more copper). All the
copper heat then has to cross, per tooth:

    conductor insulation + slot liner  ->  tooth  ->  yoke  ->  stator/jacket fit  ->  coolant film

A lumped model of that path per tooth, per metre of stack. Heat made in the end turns
runs along the copper into the slot (copper is a far better conductor than anything
else on the path), so the whole copper loss is taken through the stack length. The
end turns' extra temperature rise is estimated separately.

Every thermal number here is [ASSUMED]; the table at the end varies the ones that
matter. Torque comes from the same field solution as the hollow design (fea.py).

    python qdd-liquid/jacket.py        (or: make qdd-liquid-jacket)
"""

import sys
from dataclasses import dataclass, replace
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402
import fea  # noqa: E402
import study  # noqa: E402


@dataclass(frozen=True)
class Path_:
    k_liner: float = 0.2          # W/mK: aramid/polyimide liner and wire insulation, unpotted
    gap_mm: float = 0.0           # extra air gap between conductor and liner (0 = pressed/varnished)
    k_iron: float = 25.0          # W/mK, laminations in-plane
    h_fit: float = 2000.0         # W/m2K, stator pressed into the jacket sleeve
    h_film: float = 1500.0        # W/m2K, coolant in the 1.5 mm annular channel
    t_coolant: float = 60.0       # C


def solid(d):
    return replace(d, cond_bore=0.0)


def resistance(d, p):
    """K·m/W from the copper of one tooth's two coil sides to the coolant, and its parts."""
    r_a, r_b = (x * 1e3 for x in d.slot_radii)
    pitch = lambda r: 2 * pi * r / d.slots
    n = 2 * d.turns                                         # conductors against this tooth
    face = (d.cond_h - 2 * d.insulation) * 1e-3             # each one's face on the tooth wall
    t_liner = (d.liner + d.insulation) * 1e-3
    liner = (t_liner / p.k_liner + p.gap_mm * 1e-3 / 0.03) / (n * face)
    tooth_w = (pitch((r_a + r_b) / 2) - d.slot_w) * 1e-3
    tooth = (d.slot_depth / 2 * 1e-3) / (p.k_iron * tooth_w)
    yoke = (d.stator_yoke * 1e-3) / (p.k_iron * pitch(r_b + d.stator_yoke / 2) * 1e-3)
    arc = pitch(r_b + d.stator_yoke) * 1e-3
    fit, film = 1 / (p.h_fit * arc), 1 / (p.h_film * arc)
    parts = dict(liner=liner, tooth=tooth, yoke=yoke, fit=fit, film=film)
    return sum(parts.values()), parts


def jacket_only(d, p, t_cu_max):
    """Copper loss, current density and end-turn rise with the slot copper at t_cu_max."""
    r, _ = resistance(d, p)
    watts = (t_cu_max - p.t_coolant) / r * d.stack * 1e-3 * d.slots
    rho = dz.RHO_CU_20 * (1 + dz.ALPHA_CU * (t_cu_max - 20))
    vol = d.cu_area * d.coil_len * d.slots
    j = sqrt(watts / (rho * vol))
    # end turn: uniform heating of a bar of half the end-turn length, cooled at its root
    half = d.end_turn_mm / 2 * 1e-3
    end_rise = rho * j * j * half * half / (2 * 390.0)
    return watts, j, end_rise


def main():
    hollow, sol = dz.SHALLOW, solid(dz.SHALLOW)
    motor = fea.Motor(hollow)
    ats = np.array([fea.amp_turns_of(hollow, j) for j in study.J_LEVELS])
    torque = motor.curve(ats)
    t_of = lambda d, j: float(np.interp(fea.amp_turns_of(d, j), ats, torque))

    p = Path_()
    total, parts = resistance(sol, p)
    print(f"solid conductor: {sol.cu_area*1e6:.1f} mm² against the hollow one's {hollow.cu_area*1e6:.1f} "
          f"({sol.cu_area/hollow.cu_area*100-100:.0f}% more copper)")
    print(f"heat path per tooth, K per (W per mm of stack): " + ", ".join(f"{k} {v*1e3:.0f}" for k, v in parts.items())
          + f"  = {total*1e3:.0f}")

    print(f"\n=== jacket only, coolant at {p.t_coolant:.0f} C ===")
    print(f"  {'slot copper at':>15} {'copper loss':>12} {'A/mm²':>6} {'A rms':>6} {'torque':>8} {'end turns hotter by':>20}")
    rows = {}
    for t in (130, 155, 180, 200):
        w, j, rise = jacket_only(sol, p, t)
        rows[t] = (w, j, t_of(sol, j))
        print(f"  {t:13.0f} C {w:10.0f} W {j/1e6:6.1f} {j*sol.cu_area:6.0f} {t_of(sol, j):6.1f} N·m {rise:18.1f} K")

    print("\n=== against the hollow conductor (same iron, same slots) ===")
    print(f"  {'':34} {'torque':>8} {'heat':>7} {'hottest copper':>15}")
    for label, heat in (("hollow, 850 W", 850), ("hollow, at the 3 L/min limit", None)):
        js = np.geomspace(2e6, 140e6, 70)
        cool = [(j, dz.cooled(hollow, j)) for j in js]
        cool = [(j, c) for j, c in cool if c is not None]
        if heat:
            j = float(np.interp(heat, [dz.total_loss(c) for _, c in cool], [j for j, _ in cool]))
        else:
            j = cool[-1][0]
        c = dz.cooled(hollow, j)
        print(f"  {label:34} {t_of(hollow, j):6.1f} N·m {dz.total_loss(c):5.0f} W {c['t_max']:13.0f} C")
    for t in (130, 180):
        w, j, tq = rows[t]
        print(f"  {'jacket only, copper at ' + str(t) + ' C':34} {tq:6.1f} N·m {w:5.0f} W {t:13.0f} C")

    print("\n=== what the jacket-only figure depends on (copper at 180 C) ===")
    print(f"  {'assumption':46} {'copper loss':>12} {'torque':>8}")
    cases = [
        ("as above", p),
        ("0.05 mm air gap at the liner (loose coil)", replace(p, gap_mm=0.05)),
        ("potted, liner path at 1 W/mK", replace(p, k_liner=1.0)),
        ("poor fit to the jacket (500 W/m²K)", replace(p, h_fit=500.0)),
        ("turbulent coolant film (5000 W/m²K)", replace(p, h_film=5000.0)),
        ("potted, good film, good fit (best case)", replace(p, k_liner=1.0, h_film=5000.0, h_fit=4000.0)),
        ("coolant at 30 C", replace(p, t_coolant=30.0)),
    ]
    for label, pc in cases:
        w, j, _ = jacket_only(sol, pc, 180)
        print(f"  {label:46} {w:10.0f} W {t_of(sol, j):6.1f} N·m")


if __name__ == "__main__":
    main()
