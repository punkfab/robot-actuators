"""
Thermal field solution of the jacket-only motor, to replace jacket.py's five lumped
resistances with a solved temperature field.

One slot pitch of the shallow-slot stator in cross-section, in FEMM's heat-flow solver:
tooth and yoke iron, slot liner, the eight solid conductors of the two coil sides with
their insulation, the voids between them, and the slot opening. Copper carries a uniform
heat source; the stator's outer surface is a convection boundary (press fit and coolant
film in series); the bore and the two tooth centre-lines are adiabatic. As in jacket.py,
end-turn heat is taken to enter the slot along the copper, so the whole copper loss
crosses this section.

The problem is linear, so one solve per case gives a resistance in K per (W per metre of
stack per slot), the quantity jacket.resistance() returns. Same assumed coefficients as
jacket.py; this checks the lumping, not the coefficients.

    python qdd-liquid/thermal_check.py        (or: make qdd-liquid-thermal)
"""

import json
import sys
from dataclasses import replace
from math import asin, cos, sin, sqrt, tau
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402
import fea  # noqa: E402
import jacket as jk  # noqa: E402
import study  # noqa: E402
from femm_check import OUT, Draw, rot, run_lua, wpath  # noqa: E402

K_COPPER, K_AIR = 390.0, 0.03


def slot_grid(d, open_w):
    """Rectilinear cells of the slot interior in slot coordinates (u along the depth,
    v across), each with a material and, for copper, its conductor index."""
    sa, sb = (x * 1e3 for x in d.slot_radii)
    half, ins = d.slot_w / 2, d.insulation
    env_out = half - d.liner                                    # conductor envelope against the liner
    env_in = env_out - d.cond_w
    vs = sorted({s * x for s in (-1, 1) for x in (half, env_out, env_out - ins, env_in + ins, env_in, open_w / 2)})
    us = {sa, sa + d.liner, sb - d.liner, sb}
    rows = []
    for i in range(d.turns):
        e0 = sa + d.liner + d.cond_h * i
        rows.append((e0, e0 + d.cond_h))
        us |= {e0, e0 + ins, e0 + d.cond_h - ins, e0 + d.cond_h}
    us = sorted(round(u, 6) for u in us)
    us = [u for k, u in enumerate(us) if k == 0 or u - us[k - 1] > 1e-6]
    cells = []
    for u0, u1 in zip(us, us[1:]):
        for v0, v1 in zip(vs, vs[1:]):
            u, v = (u0 + u1) / 2, (v0 + v1) / 2
            mat, cond = "fill", None
            if abs(v) > env_out or u < sa + d.liner or u > sb - d.liner:
                mat = "liner"
            elif abs(v) > env_in and u < rows[-1][1]:
                i = next(k for k, (a, b) in enumerate(rows) if a < u < b)
                inside = env_in + ins < abs(v) < env_out - ins and rows[i][0] + ins < u < rows[i][1] - ins
                mat, cond = ("copper", (i, v > 0)) if inside else ("insulation", None)
            cells.append(dict(u=(u0, u1), v=(v0, v1), mat=mat, cond=cond))
    return us, vs, cells


def thermal_lua(d, p, k_fill, ans_out, mesh=0.25):
    m_r = dict(bore=d.slot_radii[0] * 1e3 - d.tooth_tip, outer=d.slot_radii[1] * 1e3 + d.stator_yoke)
    sa, sb = (x * 1e3 for x in d.slot_radii)
    sp_ = tau / d.slots
    rb, ro = m_r["bore"], m_r["outer"]
    open_w = 2 * rb * sin(0.2 * np.pi / d.slots)                # fea.Motor's slot opening
    us, vs, cells = slot_grid(d, open_w)
    cu_cells = [c for c in cells if c["mat"] == "copper"]
    cu_area = sum((c["u"][1] - c["u"][0]) * (c["v"][1] - c["v"][0]) for c in cu_cells) * 1e-6      # m2 per slot
    h_eff = 1 / (1 / p.h_fit + 1 / p.h_film)
    k_ins = p.k_liner if not p.gap_mm else (d.liner + d.insulation) / ((d.liner + d.insulation) / p.k_liner + p.gap_mm / K_AIR)

    g = Draw("hi")
    c = g.cmd
    c += ["newdocument(2)", 'hi_probdef("millimeters","planar",1e-8,1000,30)']
    c += [f'hi_addmaterial("iron",{p.k_iron},{p.k_iron},0,0)', f'hi_addmaterial("liner",{k_ins},{k_ins},0,0)',
          f'hi_addmaterial("insulation",{k_ins},{k_ins},0,0)', f'hi_addmaterial("fill",{k_fill},{k_fill},0,0)',
          f'hi_addmaterial("air",{K_AIR},{K_AIR},0,0)',
          f'hi_addmaterial("copper",{K_COPPER},{K_COPPER},{1.0 / cu_area:.4f},0)',     # 1 W per metre of stack per slot
          f'hi_addboundprop("jacket",2,0,0,0,{h_eff:.3f},0)']
    for u in us:
        for v0, v1 in zip(vs, vs[1:]):
            g.seg((u, v0), (u, v1))
    for v in vs:
        for u0, u1 in zip(us, us[1:]):
            g.seg((u0, v), (u1, v))
    da, ub = asin(open_w / 2 / rb), sqrt(rb * rb - open_w * open_w / 4)
    for s in (-1, 1):                                            # slot opening, bore, tooth centre-lines
        g.seg((ub, s * open_w / 2), (sa, s * open_w / 2))
        a0, a1 = sorted((s * da, s * sp_ / 2))
        g.arc(rb, a0, a1, 1.0)
        g.seg(rot(rb, 0, s * sp_ / 2), rot(ro, 0, s * sp_ / 2))
    g.arc(ro, -sp_ / 2, sp_ / 2, 0.5)
    c += [f"hi_selectarcsegment({ro:.5f},0)", 'hi_setarcsegmentprop(0.5,"jacket",0,0,"<None>")', "hi_clearselected()"]
    for cell in cells:
        g.label(sum(cell["u"]) / 2, sum(cell["v"]) / 2, cell["mat"], mesh)
    g.label((rb + sa) / 2, 0, "air", mesh)
    g.label(*rot((sb + ro) / 2, 0, sp_ * 0.4), "iron", 0.6)
    c += [f'hi_saveas("{wpath(OUT / "slot.feh")}")', "hi_analyze(1)", "hi_loadsolution()",
          f'f = openfile("{wpath(ans_out)}","w")']
    for cell in cu_cells:
        for u in np.linspace(*cell["u"], 5):
            for v in np.linspace(*cell["v"], 5):
                c.append(f'T = ho_getpointvalues({u:.5f},{v:.5f}); write(f,"C {cell["cond"][0]} {int(cell["cond"][1])} ",T,"\\n")')
    probes = dict(tooth=rot((sa + sb) / 2, 0, sp_ * 0.45), yoke=((sb + ro) / 2, 0.0), skin=(ro - 0.01, 0.0))
    for name, (x, y) in probes.items():
        c.append(f'T = ho_getpointvalues({x:.5f},{y:.5f}); write(f,"P {name} ",T,"\\n")')
    c += ["closefile(f)", "ho_close()", "quit()"]
    return "\n".join(c) + "\n", cu_area


def solve(d, p, k_fill=K_AIR, tag="slot"):
    """K per (W per metre of stack per slot): hottest copper, mean copper, and probes."""
    ans = OUT / f"{tag}.txt"
    lua, cu_area = thermal_lua(d, p, k_fill, ans)
    run_lua(lua, tag)
    cu, probes = {}, {}
    for line in ans.read_text().split("\n"):
        w = line.split()
        if w[:1] == ["C"]:
            cu.setdefault((int(w[1]), int(w[2])), []).append(float(w[3]))
        elif w[:1] == ["P"]:
            probes[w[1]] = float(w[2])
    per = {k: float(np.mean(v)) for k, v in cu.items()}
    allv = np.concatenate([np.array(v) for v in cu.values()])
    return dict(r_max=float(allv.max()), r_mean=float(np.mean(list(per.values()))), per_turn=[per[(i, 1)] for i in range(d.turns)],
                probes=probes, cu_area=cu_area)


def at_copper_temp(d, p, r, t_cu, t_of):
    """Copper loss, current density and torque with the hottest copper at t_cu, using
    resistance r in place of jacket.resistance()."""
    watts = (t_cu - p.t_coolant) / r * d.stack * 1e-3 * d.slots
    rho = dz.RHO_CU_20 * (1 + dz.ALPHA_CU * (t_cu - 20))
    j = sqrt(watts / (rho * d.cu_area * d.coil_len * d.slots))
    return watts, j, t_of(d, j)


def main():
    hollow, sol = dz.SHALLOW, jk.solid(dz.SHALLOW)
    motor = fea.Motor(hollow)
    ats = np.array([fea.amp_turns_of(hollow, j) for j in study.J_LEVELS])
    torque = motor.curve(ats)
    t_of = lambda d, j: float(np.interp(fea.amp_turns_of(d, j), ats, torque))

    p = jk.Path_()
    cases = [
        ("as modelled: liner 0.2 W/mK, voids air", p, K_AIR),
        ("loose coil, 0.05 mm air gap at the liner", replace(p, gap_mm=0.05), K_AIR),
        ("poor fit to the jacket (500 W/m²K)", replace(p, h_fit=500.0), K_AIR),
        ("turbulent coolant film (5000 W/m²K)", replace(p, h_film=5000.0), K_AIR),
        ("potted: liner path and voids at 1 W/mK", replace(p, k_liner=1.0), 1.0),
        ("potted, good film, good fit (best case)", replace(p, k_liner=1.0, h_film=5000.0, h_fit=4000.0), 1.0),
    ]
    out = []
    print(f"  {'case (copper at 180 C, coolant 60 C)':42} {'lumped':>18} {'field, hottest copper':>24} {'field, mean copper':>20}")
    for k, (label, pc, k_fill) in enumerate(cases):
        r_l, parts = jk.resistance(sol, pc)
        f = solve(sol, pc, k_fill, tag=f"slot{k}")
        w_l, _, tq_l = at_copper_temp(sol, pc, r_l, 180, t_of)
        w_f, _, tq_f = at_copper_temp(sol, pc, f["r_max"], 180, t_of)
        w_m, _, tq_m = at_copper_temp(sol, pc, f["r_mean"], 180, t_of)
        out.append(dict(label=label, r_lumped=r_l, parts=parts, **f, w_lumped=w_l, t_lumped=tq_l, w_field=w_f, t_field=tq_f,
                        w_field_mean=w_m, t_field_mean=tq_m))
        print(f"  {label:42} {tq_l:6.1f} N·m {w_l:5.0f} W {tq_f:10.1f} N·m {w_f:5.0f} W {tq_m:8.1f} N·m {w_m:5.0f} W")
    b = out[0]
    print(f"\nresistance, K per (W per mm of stack per slot): lumped {b['r_lumped']*1e3:.0f}, "
          f"field {b['r_max']*1e3:.0f} to the hottest copper, {b['r_mean']*1e3:.0f} to the mean")
    print("as modelled, rise above coolant per turn (bore side first), share of the hottest: "
          + ", ".join(f"{x/b['r_max']*100:.0f}%" for x in b["per_turn"]))
    print("as modelled, share of the rise reached at: " + ", ".join(f"{k} {v/b['r_max']*100:.0f}%" for k, v in b["probes"].items()))
    (OUT / "thermal.json").write_text(json.dumps(out, indent=1))


if __name__ == "__main__":
    main()
