"""
Figures and numbers for the paper (paper/main.tex), from the same code as the README.

  figures/torque.pdf   shallow-slot design: torque against current density, and torque
                       against copper loss for hollow conductors and for the jacket alone
  figures/field.png    iron flux density at three currents
  numbers.json         every figure the paper quotes that is not already in the README,
                       including the two solver checks (mesh refinement, low-current
                       torque against the closed form)

    python qdd-liquid/paper/figures.py        (or: make qdd-liquid-paper; a few minutes)
"""

import json
import sys
from math import sqrt
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import design as dz  # noqa: E402
import fea  # noqa: E402
import jacket as jk  # noqa: E402
import study  # noqa: E402

FIG = HERE / "figures"
D = dz.SHALLOW


def mesh_check(j=32.5e6):
    """Torque at one current and angle on three meshes."""
    at = fea.amp_turns_of(D, j)
    rows, g = [], None
    for h_gap, h_coarse in ((0.35, 2.2), (0.25, 1.6), (0.15, 1.0)):
        m = fea.Motor(D, h_gap=h_gap, h_coarse=h_coarse)
        if g is None:
            g, _ = m.best_gamma(at)
        rows.append(dict(h_gap=h_gap, h_coarse=h_coarse, triangles=int(len(m.tri)),
                         torque=float(m.solve(at, g)["torque"]), b1=float(m.b1_gap())))
    return rows


def closed_form_check(motor, j=3e6):
    """Low-current torque against 3 kw N B1 r L I with the solver's own air-gap field."""
    at = fea.amp_turns_of(D, j)
    g, _ = motor.best_gamma(at)
    b1 = motor.b1_gap()
    r_mid = sum(motor.gap) / 2 * 1e-3
    formula = 3 * D.kw * (D.slots // 3) * b1 * r_mid * D.stack * 1e-3 * at
    return dict(j=j, b1=float(b1), fea=float(motor.solve(at, g)["torque"]), formula=float(formula))


def main():
    FIG.mkdir(exist_ok=True)
    v = study.Variant("shallow", D, keep_turns=True)
    sol = jk.solid(D)
    p = jk.Path_()
    t_sol = lambda j: float(np.interp(fea.amp_turns_of(sol, j), [fea.amp_turns_of(D, x) for x in study.J_LEVELS], v.torque))

    # hollow: copper loss against torque up to the cooling limit
    cool = [dz.cooled(D, j) for j in v.j_ok]
    cu_h = np.array([c["p_cu"] for c in cool])
    tq_h = np.array([v.t_at_j(j) for j in v.j_ok])
    # jacket only: slot copper from 80 to 200 C
    temps = np.arange(80, 201, 10)
    # resistance from the thermal field solution (thermal_check.py) in place of the lumped one
    r_field = json.loads((HERE.parent / "out" / "femm" / "thermal.json").read_text())[0]["r_max"]
    scale = jk.resistance(sol, p)[0] / r_field
    jrows = [(w * scale, j * sqrt(scale), rise * scale) for w, j, rise in (jk.jacket_only(sol, p, t) for t in temps)]
    cu_j = np.array([r[0] for r in jrows])
    tq_j = np.array([t_sol(r[1]) for r in jrows])
    air = dz.air_cooled(D, 88.0)                 # current from the sizing model, torque from the field solution
    air["torque"] = float(np.interp(sqrt(2) * air["j"] * D.slot_area * 0.45 / 2,
                                    [fea.amp_turns_of(D, x) for x in study.J_LEVELS], v.torque))

    jd = v.j_kt_drop()
    c850 = dz.cooled(D, v.at_heat(850)[1])
    drv = lambda i: 3 * i * i * (D.r_dson_mohm + D.r_lead_mohm) * 1e-3 + 3 * D.v_bus * 0.9 * i * D.t_sw * D.f_sw
    numbers = dict(
        mesh=mesh_check(), closed_form=closed_form_check(v.motor),
        curve=[dict(j=float(j), torque=float(t), kt=float(v.kt_ratio(j))) for j, t in zip(study.J_LEVELS, v.torque)],
        j_cool=float(v.j_cool), t_cool=float(v.t_at_j(v.j_cool)), j_kt20=jd, t_kt20=v.t_at_j(jd) if jd else None,
        at_850=dict(torque=v.at_heat(850)[0], j=v.at_heat(850)[1], **{k: float(c850[k]) for k in
                    ("i_rms", "p_cu", "p_drv", "p_pump", "p_fe", "t_max", "t_out", "t_skin", "t_mag", "flow_lpm", "dp_bar", "re", "v_ph")}),
        jacket=[dict(t_cu=int(t), p_cu=float(r[0]), j=float(r[1]), i_rms=float(r[1] * sol.cu_area), torque=float(q),
                     p_drv=float(drv(r[1] * sol.cu_area)), end_rise=float(r[2]))
                for t, r, q in zip(temps, jrows, tq_j)],
        air=dict(torque=float(air["torque"]), j=float(air["j"])),
        cu_area_hollow=D.cu_area * 1e6, cu_area_solid=sol.cu_area * 1e6, mass=D.mass(), bore_free=D.bore_free_mm,
        r_gap=D.r_gap * 1e3, n_phase=D.n_phase, triangles=int(len(v.motor.tri)),
    )
    (HERE / "numbers.json").write_text(json.dumps(numbers, indent=1))

    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.4, 3.1))
    jm = study.J_LEVELS / 1e6
    a.plot(jm, v.kt0 * study.J_LEVELS, "--", color="0.55", lw=1, label="no saturation")
    a.plot(jm, v.torque, "o-", color="#1f4e79", ms=3.5, lw=1.6, label="field solution")
    a.axvline(v.j_cool / 1e6, color="#2a7f3f", lw=1)
    a.text(v.j_cool / 1e6 + 2, 2, "cooling limit\n(3 L/min)", color="#2a7f3f", fontsize=7.5)
    if jd:
        a.plot([jd / 1e6], [v.t_at_j(jd)], "s", color="#b03a2e", ms=5)
        a.annotate("torque per amp\ndown 20%", (jd / 1e6, v.t_at_j(jd)), (jd / 1e6 + 8, v.t_at_j(jd) - 9),
                   fontsize=7.5, color="#b03a2e", arrowprops=dict(arrowstyle="-", color="#b03a2e", lw=0.7))
    a.set(xlabel="copper current density (A/mm$^2$ rms)", ylabel="motor torque (N·m)", xlim=(0, 145), ylim=(0, 45))
    a.legend(frameon=False, fontsize=7.5, loc="upper left")

    b.plot(cu_h, tq_h, color="#1f4e79", lw=1.6, label="hollow conductors")
    b.plot(cu_j, tq_j, color="#c87a2f", lw=1.6, label="solid, jacket only")
    for t in (130, 180):
        k = list(temps).index(t)
        b.plot([cu_j[k]], [tq_j[k]], "o", color="#c87a2f", ms=4)
        b.annotate(f"copper {t} °C", (cu_j[k], tq_j[k]), (cu_j[k] + 70, tq_j[k] - (3.6 if t == 130 else 2.2)), fontsize=7.5, color="#c87a2f")
    k = int(np.argmin(np.abs(cu_h - c850["p_cu"])))
    b.plot([cu_h[k]], [tq_h[k]], "o", color="#1f4e79", ms=4)
    b.annotate(f"copper {c850['t_max']:.0f} °C", (cu_h[k], tq_h[k]), (cu_h[k] + 60, tq_h[k] - 2.6), fontsize=7.5, color="#1f4e79")
    b.plot([88], [air["torque"]], "x", color="0.3", ms=6)
    b.annotate("air-cooled, 88 W", (88, air["torque"]), (170, air["torque"] - 2.6), fontsize=7.5, color="0.3")
    b.set(xlabel="copper loss (W)", ylabel="motor torque (N·m)", xlim=(0, 1800), ylim=(0, 24))
    b.legend(frameon=False, fontsize=7.5, loc="lower right")
    fig.tight_layout()
    fig.savefig(FIG / "torque.pdf")

    g, _ = v.motor.best_gamma(fea.amp_turns_of(D, 3e6))
    levels = (0.0, numbers["at_850"]["j"], jd or 67e6)
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.5))
    for ax, j in zip(axes, levels):
        pc, _ = study.field_plot(ax, v.motor, fea.amp_turns_of(D, j), g, f"{j/1e6:.0f} A/mm$^2$")
        ax.title.set_fontsize(9)
    fig.colorbar(pc, ax=list(axes), shrink=0.8, label="|B| in iron (T)", pad=0.02)
    fig.savefig(FIG / "field.png", bbox_inches="tight", dpi=300)
    print(json.dumps({k: numbers[k] for k in ("mesh", "closed_form", "j_cool", "t_cool", "j_kt20", "t_kt20", "at_850", "air")}, indent=1))
    for r in numbers["jacket"]:
        print(r)


if __name__ == "__main__":
    main()
