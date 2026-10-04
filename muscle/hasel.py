"""
A HASEL cell beside the SMA one: one Peano-HASEL pouch, stacked into a muscle unit.

Option A from README.md: the cell technology of Bettini et al. 2026 (Nat. Commun.,
doi 10.1038/s41467-026-77664-0), so their control recipe can be tried where it was
meant to work and compared with our SMA cell on the same task.

Geometry and force are the quasi-static pouch model of Kellaris et al. 2019 (Extreme
Mech. Lett. 29, 100449), as used in the paper's Eq. 11–12 and Supplementary Eq. 9–12:

    a pouch of film length Lp holds a fixed cross-section A of oil. Electrodes zip
    shut over a length le from one edge; the rest, lp = Lp - le, bulges into two
    circular arcs of half-angle alpha. Zipping forces alpha up, which shortens the pouch.

        lp(alpha) = sqrt(2 A alpha^2 / (alpha - sin alpha cos alpha))
        length    = le + lp sin(alpha) / alpha
        C         = eps0 eps_r w le / (2 t)
        F         = 1/2 V^2 dC/dx = (eps0 eps_r w V^2 / 4t) * cos(alpha) / (1 - cos(alpha))

Sources, tagged in the spec:
  [PAPER]    stated in Bettini et al.: BoPET shell, 15 pouches in series, up to 8 kV
  [KELLARIS] the model above; the oil fill that makes the pouch a cylinder when the
             electrodes are fully zipped
  [ASSUMED]  everything dimensional (film thickness, pouch and electrode length,
             width) and every dynamic number. Neither the paper nor its supplement
             gives them. Dimensions are typical of published Peano-HASELs.

Dynamics (HaselPlant): the paper's Eq. 13 in its overdamped limit, c x' = F_es - load,
with charging treated as instant. c is set from an assumed time constant.

    python muscle/hasel.py        (or: make muscle-hasel)
"""

import sys
from dataclasses import dataclass, replace
from math import ceil, cos, pi, sin, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from cell import SMACellSpec  # noqa: E402

EPS0 = 8.854e-12


@dataclass
class HaselCellSpec:
    """One Peano-HASEL pouch, and how many are stacked."""

    pouch_length_mm: float = 20.0       # Lp, along the pull direction        [ASSUMED]
    electrode_fraction: float = 0.5     # Le / Lp                              [ASSUMED]
    width_mm: float = 50.0              # w, across the pull direction         [ASSUMED]
    film_um: float = 18.0               # t, one film                          [ASSUMED]
    eps_r: float = 3.3                  # BoPET                                [PAPER: material]
    v_max: float = 8000.0               # V                                    [PAPER]
    series_count: int = 15              # pouches in series -> stroke x        [PAPER]
    parallel_count: int = 1             # stacks side by side -> force x
    tau_s: float = 0.05                 # mechanical time constant, mid-stroke [ASSUMED]

    # ---- geometry ----------------------------------------------------------
    @property
    def lp0(self) -> float:
        return self.pouch_length_mm * 1e-3

    @property
    def le_max(self) -> float:
        return self.lp0 * self.electrode_fraction

    @property
    def area(self) -> float:
        """Oil cross-section: a full cylinder when the electrodes are zipped shut. [KELLARIS]"""
        return (self.lp0 - self.le_max) ** 2 / pi

    def lp(self, alpha):
        return np.sqrt(2 * self.area * alpha ** 2 / (alpha - np.sin(alpha) * np.cos(alpha)))

    @property
    def alpha0(self) -> float:
        """Rest half-angle: lp(alpha0) = Lp, by bisection (lp falls as alpha rises)."""
        lo, hi = 1e-3, pi / 2
        for _ in range(60):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self.lp(mid) > self.lp0 else (lo, mid)
        return 0.5 * (lo + hi)

    def stroke(self, alpha):
        """Contraction of one pouch from rest, m."""
        a0 = self.alpha0
        lp = self.lp(alpha)
        return self.lp0 * sin(a0) / a0 - (lp * np.sin(alpha) / alpha + self.lp0 - lp)

    def capacitance(self, alpha):
        return EPS0 * self.eps_r * self.width_mm * 1e-3 * (self.lp0 - self.lp(alpha)) / (2 * self.film_um * 1e-6)

    def force(self, alpha, v):
        """Contractile force of one pouch at half-angle alpha and voltage v, N."""
        k = EPS0 * self.eps_r * self.width_mm * 1e-3 * v ** 2 / (4 * self.film_um * 1e-6)
        return k * np.cos(alpha) / (1 - np.cos(alpha))

    def table(self, n=4000):
        """(stroke, force per V^2) over the whole travel, stroke ascending."""
        alpha = np.linspace(self.alpha0, pi / 2, n)
        return self.stroke(alpha), self.force(alpha, 1.0)

    # ---- the stacked unit ----------------------------------------------------
    @property
    def free_stroke_mm(self) -> float:
        return float(self.stroke(pi / 2)) * 1e3 * self.series_count

    @property
    def free_strain(self) -> float:
        return float(self.stroke(pi / 2)) / (self.lp0 * sin(self.alpha0) / self.alpha0)

    @property
    def length_mm(self) -> float:
        return self.pouch_length_mm * self.series_count

    def unit_force(self, stroke_mm, v=None):
        """Force of the whole unit at a total stroke, N."""
        x, f = self.table()
        per = stroke_mm * 1e-3 / self.series_count
        return float(np.interp(per, x, f)) * (self.v_max if v is None else v) ** 2 * self.parallel_count

    @property
    def field_v_per_um(self) -> float:
        return self.v_max / (2 * self.film_um)

    @property
    def energy_per_stroke_j(self) -> float:
        """1/2 C V^2 with every electrode fully zipped: the most one stroke can cost."""
        return 0.5 * float(self.capacitance(pi / 2)) * self.v_max ** 2 * self.series_count * self.parallel_count

    def size_for(self, target_force_n, target_stroke_mm) -> "HaselCellSpec":
        """Copy sized to give the force AT the stroke, at v_max: series from stroke
        (each pouch used to at most a third of its free stroke), then parallel stacks."""
        per_max = float(self.stroke(pi / 2)) * 1e3 / 3
        n_ser = max(self.series_count, ceil(target_stroke_mm / per_max))
        one = replace(self, series_count=n_ser, parallel_count=1)
        return replace(one, parallel_count=max(1, ceil(target_force_n / one.unit_force(target_stroke_mm))))

    def report(self):
        print(f"\n=== HASEL pouch (BoPET {self.film_um:.0f} um, {self.pouch_length_mm:.0f} x "
              f"{self.width_mm:.0f} mm, electrodes {self.electrode_fraction*100:.0f}%) ===")
        print(f"rest half-angle ... {np.degrees(self.alpha0):.1f} deg; free strain {self.free_strain*100:.1f}%")
        print(f"field at {self.v_max/1e3:.0f} kV ..... {self.field_v_per_um:.0f} V/um across the two films")
        print(f"blocked force ..... {float(self.force(self.alpha0, self.v_max)):.1f} N per pouch "
              f"(falls to 0 at full zip)")
        print(f"\nunit: {self.parallel_count} parallel x {self.series_count} series, "
              f"{self.length_mm:.0f} mm long")
        print(f"  free stroke ..... {self.free_stroke_mm:.1f} mm")
        for s in (0.0, 0.25, 0.5):
            mm = s * self.free_stroke_mm
            print(f"  force at {mm:5.1f} mm  {self.unit_force(mm):6.1f} N")
        print(f"  energy/stroke ... {self.energy_per_stroke_j*1e3:.0f} mJ at most (1/2 C V^2)")


class HaselPlant:
    """The unit pulling a tendon against a preload and a spring; input u = V / v_max.

    Overdamped: c x' = F_es(x, V_eff) - (f0 + k x). c comes from spec.tau_s at
    mid-stroke and half the squared voltage. Output is stroke over the stroke reached
    at v_max. `loop` is a play operator on u (half-width, fraction of v_max) standing
    in for the hysteresis the paper reports but doesn't quantify. [ASSUMED]"""

    def __init__(self, spec, f0, f1, stroke_mm, loop=0.0, dt=0.001):
        self.spec, self.loop, self.dt, self.tau = spec, loop, dt, spec.tau_s
        xs, fs = spec.table()
        self.x_full = stroke_mm * 1e-3
        self.n = 2000
        self.dx = self.x_full * 1.2 / self.n                 # uniform grid for a fast lookup
        grid = np.arange(self.n + 2) * self.dx
        self.ftab = (np.interp(grid / spec.series_count, xs, fs)
                     * spec.v_max ** 2 * spec.parallel_count).tolist()
        self.f0, self.k = f0, (f1 - f0) / self.x_full
        i = int(0.5 * self.x_full / self.dx)
        stiff = self.k - 0.5 * (self.ftab[i + 1] - self.ftab[i - 1]) / (2 * self.dx)
        self.c = self.tau * stiff
        self.x_end = self.settle(1.0)

    def f_es(self, x, u):
        i = min(max(x / self.dx, 0.0), self.n)
        j = int(i)
        return (self.ftab[j] + (i - j) * (self.ftab[j + 1] - self.ftab[j])) * u * u

    def settle(self, u):
        x = 0.0
        for _ in range(int(40 * self.tau / self.dt)):
            x = max(0.0, x + (self.f_es(x, u) - self.f0 - self.k * x) / self.c * self.dt)
        return x

    def u_range(self):
        """From just under the voltage where zipping starts, to full voltage."""
        return np.array([0.9 * sqrt(self.f0 / self.ftab[0]), 1.0])

    def start(self, a0):
        """State at output a0, approached from below."""
        lo, hi = 0.0, 1.0
        for _ in range(30):
            mid = 0.5 * (lo + hi)
            lo, hi = (mid, hi) if self.settle(mid) < a0 * self.x_end else (lo, mid)
        return a0 * self.x_end, hi - self.loop

    def rest(self):
        return 0.0, 0.0

    def act(self, state):
        return state[0] / self.x_end

    def step(self, state, u):
        x, ue = state
        ue = min(max(ue, u - self.loop), u + self.loop)
        x = max(0.0, x + (self.f_es(x, max(ue, 0.0)) - self.f0 - self.k * x) / self.c * self.dt)
        return x, ue


UNIT_FORCE_N, UNIT_STROKE_MM = 10.0, 10.0       # the 'muscle unit' cell.py sizes for


def unit():
    return HaselCellSpec().size_for(UNIT_FORCE_N, UNIT_STROKE_MM)


def plant(loop=0.0, tau=None, dt=0.001):
    """The sized unit against a load rising from half to full unit force over the stroke."""
    spec = unit() if tau is None else replace(unit(), tau_s=tau)
    return HaselPlant(spec, UNIT_FORCE_N / 2, UNIT_FORCE_N, UNIT_STROKE_MM, loop, dt)


if __name__ == "__main__":
    base = HaselCellSpec()
    # the closed-form force against 1/2 V^2 dC/dx taken numerically
    a = np.linspace(base.alpha0 * 1.05, 1.5, 9)
    h = 1e-6
    num = 0.5 * base.v_max ** 2 * (base.capacitance(a + h) - base.capacitance(a - h)) \
        / (base.stroke(a + h) - base.stroke(a - h))
    err = np.max(np.abs(num / base.force(a, base.v_max) - 1))
    print(f"force formula vs 1/2 V^2 dC/dx: max difference {err*100:.3f}%")

    h_unit = unit()
    h_unit.report()

    sma = SMACellSpec().size_for(UNIT_FORCE_N, UNIT_STROKE_MM)
    p = plant()
    print(f"\n=== the same {UNIT_FORCE_N:.0f} N / {UNIT_STROKE_MM:.0f} mm unit, both cells ===")
    print(f"  {'':22} {'SMA':>14} {'HASEL':>14}")
    rows = [
        ("cells", f"{sma.parallel_count} x {sma.series_count} wires",
         f"{h_unit.parallel_count} x {h_unit.series_count} pouches"),
        ("active length", f"{sma.active_length_mm * sma.series_count:.0f} mm", f"{h_unit.length_mm:.0f} mm"),
        ("strain used", f"{sma.usable_strain*100:.0f}%", f"{UNIT_STROKE_MM / h_unit.length_mm * 100:.1f}%"),
        ("drive", f"{sma.rec_current_a:.2f} A, low V", f"{h_unit.v_max/1e3:.0f} kV, ~0 A"),
        ("time constant", f"{sma.tau_cool_s:.1f} s (cooling)", f"{h_unit.tau_s*1e3:.0f} ms [ASSUMED]"),
        ("hold power", f"{sma.hold_power_w:.1f} W", "~0 (leakage not modelled)"),
        ("energy per stroke", f"{sma.hold_power_w * sma.tau_heat_s:.1f} J",
         f"{h_unit.energy_per_stroke_j*1e3:.0f} mJ max"),
    ]
    for name, s, hh in rows:
        print(f"  {name:22} {s:>14} {hh:>14}")
    print(f"\n  HASEL zipping starts at {p.u_range()[0] / 0.9 * h_unit.v_max / 1e3:.1f} kV against the "
          f"{p.f0:.0f} N preload; reaches {p.x_end*1e3:.1f} mm at {h_unit.v_max/1e3:.0f} kV")
