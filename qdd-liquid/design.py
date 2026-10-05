"""
A direct-liquid-cooled, high-current QDD motor: first-order sizing.

Thesis under test: solid-state drives now handle very high current efficiently, so
the limit on a robot joint motor is heat. Take the heat out actively, run the copper
hard, and you need less gear reduction.

The design this file sizes:
  - hollow rectangular copper conductors, coolant flowing through the bore of every turn
  - few turns of large conductor: a low-voltage, very-high-current machine
  - the drive on the motor, on the same coolant loop
  - an outer jacket fed with the coldest coolant, so the skin stays near loop-return
    temperature whatever the copper is doing
  - flow order ("regenerative", as in a rocket nozzle: the cold fluid lines the
    outside before it reaches the heat source):
        radiator -> drive cold plate -> outer jacket -> conductor bores -> radiator

What the model does, per design:
  electromagnetics   T = 3 kw N_ph B1 r_g L I_peak, with a saturation knee [ASSUMED]
  copper loss        rho(T) J^2 V, end turns included
  hydraulics         bore velocity from the pump pressure budget (laminar/turbulent)
  heat transfer      bulk rise + film rise along each coil's bore
  solve              the current density at which the hottest copper hits its limit
  drive              conduction (R_ds(on)), leads, switching

The saturation knee and air-gap field here are placeholders: fea.py solves the field
and study.py uses its torque instead (the knee turned out pessimistic).

What it does not do: field solution (that is fea.py),
demagnetization, AC loss beyond a skin-depth check, coil-bend effects on flow,
manifold losses, mechanical design. Everything [ASSUMED] is a number to replace.

Baseline for comparison: RobStride RS04 (robstride/specs.py), the same 120 mm class.

    python qdd-liquid/design.py        (or: make qdd-liquid)
"""

import sys
from dataclasses import dataclass, replace
from math import exp, log, pi, sqrt
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE.parent / "robstride"))
from specs import SPECS  # noqa: E402

RHO_CU_20 = 1.72e-8        # ohm m
ALPHA_CU = 0.00393         # 1/K
MU0 = 4e-7 * pi


# ---- coolants ------------------------------------------------------------------
# (T C, density kg/m3, cp J/kgK, conductivity W/mK, viscosity Pa s). [APPROX] handbook
# values; viscosity is interpolated in log space.
COOLANTS = {
    # dielectric: can touch live copper, so coils can be electrically in series and
    # hydraulically in parallel with no insulating breaks
    "oil": dict(dielectric=True, label="dielectric oil (PAO class)", pts=[
        (40, 800, 2100, 0.140, 4.0e-3), (100, 765, 2300, 0.135, 1.3e-3), (150, 735, 2500, 0.130, 0.66e-3)]),
    # conductive: needs deionized loop or insulating breaks between coils
    "glycol": dict(dielectric=False, label="water-glycol 50/50", pts=[
        (20, 1071, 3280, 0.38, 3.9e-3), (60, 1048, 3450, 0.40, 1.4e-3), (90, 1027, 3570, 0.41, 0.8e-3)]),
}


def props(name, t):
    pts = COOLANTS[name]["pts"]
    ts = [p[0] for p in pts]
    t = min(max(t, ts[0]), ts[-1])
    f = lambda i: float(np.interp(t, ts, [p[i] for p in pts]))
    mu = exp(float(np.interp(t, ts, [log(p[4]) for p in pts])))
    return dict(rho=f(1), cp=f(2), k=f(3), mu=mu)


@dataclass(frozen=True)
class Design:
    # ---- envelope and radial build, mm ----------------------------------------
    od: float = 120.0                 # RS04's diameter
    stack: float = 25.0               # active length
    outrunner: bool = False           # False: stator outside, bonded to the jacket
    jacket: float = 3.5               # wall + channel + wall
    stator_yoke: float = 5.0
    slot_depth: float = 14.0
    tooth_tip: float = 1.5
    air_gap: float = 0.7
    magnet: float = 3.5
    rotor_yoke: float = 4.0
    rotor_shell: float = 2.5          # outrunner only: spinning can under the jacket

    # ---- electromagnetics -----------------------------------------------------
    slots: int = 24
    poles: int = 22
    kw: float = 0.95                  # winding factor, 24s/22p concentrated
    slot_frac: float = 0.5            # slot width / slot pitch (wedge slots)
    slot_w: float = 0.0               # mm; > 0 makes the slots parallel-sided, this wide
    b1_20c: float = 0.90              # air-gap fundamental, peak T, magnets at 20 C [ASSUMED]
    br_tempco: float = -0.0011        # 1/K, NdFeB
    sigma_sat_kpa: float = 90.0       # air-gap shear where teeth saturate [ASSUMED]
    kt_drop: float = 0.20             # rated point: torque per amp down this much [ASSUMED]

    # ---- winding: hollow rectangular conductor --------------------------------
    turns: int = 3                    # per coil; 2 coil sides share each slot
    pack: float = 0.70                # conductor envelopes / slot area (liner, insulation, corners)
    bore_frac: float = 0.25           # bore area / conductor envelope
    end_turn_mm: float = 18.0         # one end, per turn
    # A specific conductor instead of pack/bore_frac: insulated pitch across and along
    # the slot, and its bore. Round = tube of diameter cond_w.
    cond_w: float = 0.0               # mm, across the slot (0 = use pack and bore_frac)
    cond_h: float = 0.0               # mm, along the slot depth
    cond_bore: float = 0.0            # mm
    cond_round: bool = False
    insulation: float = 0.05          # mm per side
    liner: float = 0.25               # mm, slot liner
    t_cu_max: float = 200.0           # hottest copper, C (class 220/240 insulation, margin)

    # ---- coolant loop ---------------------------------------------------------
    coolant: str = "oil"
    t_in: float = 60.0                # C, from the radiator
    dp_bar: float = 3.0               # most the pump can put across a coil bore
    pump_eff: float = 0.4
    max_flow_lpm: float = 3.0         # what a joint-sized pump and hoses carry [ASSUMED]

    # ---- drive, on the motor --------------------------------------------------
    v_bus: float = 48.0
    r_dson_mohm: float = 0.5          # per switch position (paralleled FETs) [ASSUMED]
    r_lead_mohm: float = 0.2          # per phase, coil ends to the half-bridge [ASSUMED]
    f_sw: float = 40e3
    t_sw: float = 20e-9               # effective switching transition [ASSUMED]

    # ---- operating point for speed-dependent terms -----------------------------
    rpm: float = 300.0                # motor speed

    # ---- geometry -------------------------------------------------------------
    @property
    def r_gap(self):
        r = self.od / 2 - self.jacket
        if self.outrunner:
            return (r - self.rotor_shell - self.rotor_yoke - self.magnet - self.air_gap) * 1e-3
        return (r - self.stator_yoke - self.slot_depth - self.tooth_tip - self.air_gap) * 1e-3

    @property
    def slot_radii(self):
        """(inner, outer) radius of the slot annulus, m."""
        rg = self.r_gap * 1e3
        if self.outrunner:
            top = rg - self.air_gap - self.tooth_tip
            return (top - self.slot_depth) * 1e-3, top * 1e-3
        bot = rg + self.air_gap + self.tooth_tip
        return bot * 1e-3, (bot + self.slot_depth) * 1e-3

    @property
    def bore_free_mm(self):
        """Diameter left empty at the centre (room for a gear stage)."""
        if self.outrunner:
            return 2 * (self.slot_radii[0] * 1e3 - self.stator_yoke)
        return 2 * (self.r_gap * 1e3 - self.magnet - self.rotor_yoke)

    @property
    def slot_area(self):
        if self.slot_w:
            return self.slot_w * self.slot_depth * 1e-6
        ri, ro = self.slot_radii
        return self.slot_frac * pi * (ro ** 2 - ri ** 2) / self.slots

    @property
    def envelope_area(self):
        if self.cond_w:
            return self.cond_w * self.cond_h * 1e-6
        return self.slot_area * self.pack / (2 * self.turns)

    @property
    def bore_d(self):
        if self.cond_w:
            return self.cond_bore * 1e-3
        return sqrt(4 * self.bore_frac * self.envelope_area / pi)

    @property
    def cu_area(self):
        if not self.cond_w:
            return self.envelope_area * (1 - self.bore_frac)
        w, h = self.cond_w - 2 * self.insulation, self.cond_h - 2 * self.insulation
        bore = pi * self.cond_bore ** 2 / 4
        if self.cond_round:
            return (pi * w * w / 4 - bore) * 1e-6
        return (w * h - bore - 0.077) * 1e-6               # 0.3 mm corner radii

    def coil_grid(self):
        """(columns, rows) of this conductor that fit one coil side: half the slot
        width less liner, the slot depth less liner. Needs slot_w and cond_w."""
        cols = int((self.slot_w / 2 - self.liner - 0.05) // self.cond_w)
        rows = int((self.slot_depth - 2 * self.liner) // self.cond_h)
        return cols, rows

    @property
    def turn_len(self):
        return 2 * (self.stack + self.end_turn_mm) * 1e-3

    @property
    def coil_len(self):
        return self.turns * self.turn_len

    @property
    def n_phase(self):
        return self.slots // 3 * self.turns          # all coils of a phase in series

    def b1(self, t_mag):
        return self.b1_20c * (1 + self.br_tempco * (t_mag - 20))

    def kt_peak(self, t_mag=20.0):
        """Unsaturated torque per peak phase amp."""
        return 3 * self.kw * self.n_phase * self.b1(t_mag) * self.r_gap * self.stack * 1e-3

    def mass(self):
        """Motor only, kg: active parts x 1.3 for bearings, ends and fasteners. [ASSUMED]"""
        L = self.stack * 1e-3
        ri, ro = self.slot_radii
        cu = self.cu_area * self.coil_len * self.slots * 8960
        teeth = (pi * (ro ** 2 - ri ** 2) - self.slots * self.slot_area) * L * 7650
        if self.outrunner:
            yoke = pi * (ri ** 2 - (ri - self.stator_yoke * 1e-3) ** 2) * L * 7650
        else:
            yoke = pi * ((ro + self.stator_yoke * 1e-3) ** 2 - ro ** 2) * L * 7650
        rg = self.r_gap
        mag = 2 * pi * rg * self.magnet * 1e-3 * L * 0.85 * 7500
        ryoke = 2 * pi * rg * self.rotor_yoke * 1e-3 * L * 7650
        shell = pi * self.od * 1e-3 * (L + 2 * self.end_turn_mm * 1e-3) * self.jacket * 1e-3 * 0.6 * 2700
        return 1.3 * (cu + teeth + yoke + mag + ryoke + shell)


# The shallow-slot design (study.py: shallower slots let the rotor grow), with a
# specific conductor: parallel-sided slots 5.8 x 10 mm, each coil one column of four
# turns of hollow rectangular copper (2.5 x 2.3 mm pitch, Ø1.2 mm bore). The pitch across
# the slot leaves 0.4 mm between the two coils sharing it: cad.py found 2.6 mm clashing
# where the end bends leave the slot.
SHALLOW = Design(slot_depth=10.0, slot_w=5.8, turns=4, cond_w=2.5, cond_h=2.3, cond_bore=1.2,
                 end_turn_mm=15.0)


# ---- hydraulics and heat transfer in one coil bore -------------------------------
def friction(re):
    lam, turb = 64 / max(re, 1e-9), 0.316 * max(re, 1.0) ** -0.25
    w = min(max((re - 2300) / 1700, 0.0), 1.0)
    return (1 - w) * lam + w * turb


def nusselt(re, pr):
    turb = 0.023 * max(re, 1.0) ** 0.8 * pr ** 0.4
    w = min(max((re - 2300) / 1700, 0.0), 1.0)
    return (1 - w) * 4.36 + w * turb


def bore_flow(d, t_fluid, dp_bar):
    """Velocity that dp_bar drives through one coil. Returns v, Re, h, properties."""
    p = props(d.coolant, t_fluid)
    dia, ell = d.bore_d, d.coil_len
    lo, hi = 1e-4, 60.0
    for _ in range(60):
        v = 0.5 * (lo + hi)
        re = p["rho"] * v * dia / p["mu"]
        dp = friction(re) * ell / dia * 0.5 * p["rho"] * v * v
        lo, hi = (v, hi) if dp < dp_bar * 1e5 else (lo, v)
    re = p["rho"] * v * dia / p["mu"]
    pr = p["cp"] * p["mu"] / p["k"]
    return v, re, nusselt(re, pr) * p["k"] / dia, p


def evaluate(d, j, dp_bar=None):
    """Steady state at copper current density j (A/m2, rms) and pump pressure."""
    dp_bar = d.dp_bar if dp_bar is None else dp_bar
    i_rms = j * d.cu_area
    vol_cu = d.cu_area * d.coil_len * d.slots
    f_e = d.poles / 2 * d.rpm / 60
    iron_kg = d.mass() * 0.35
    p_fe = 2.5 * (f_e / 50) ** 1.5 * iron_kg                       # ~M19 at 1 T [APPROX]
    p_drv = 3 * i_rms ** 2 * (d.r_dson_mohm + d.r_lead_mohm) * 1e-3 \
        + 3 * d.v_bus * 0.9 * i_rms * d.t_sw * d.f_sw
    t_mean, t_max, t_fluid = d.t_in + 40, d.t_in + 60, d.t_in + 20
    for _ in range(40):
        rho_cu = RHO_CU_20 * (1 + ALPHA_CU * (t_mean - 20))
        p_cu = rho_cu * j * j * vol_cu
        v, re, h, p = bore_flow(d, t_fluid, dp_bar)
        mdot = p["rho"] * v * pi * d.bore_d ** 2 / 4 * d.slots       # coils in parallel
        t_coil_in = d.t_in + (p_drv + p_fe) / (mdot * p["cp"])       # drive plate + jacket first
        rise = p_cu / (mdot * p["cp"])
        film = (p_cu / d.slots) / (h * pi * d.bore_d * d.coil_len)
        t_fluid = t_coil_in + rise / 2
        t_mean = 0.5 * t_mean + 0.5 * (t_fluid + film)
        t_max = t_coil_in + rise + film
    t_mag = d.t_in + 0.5 * (t_mean - d.t_in)                         # between jacket and copper [ASSUMED]
    t_lin = d.kt_peak(t_mag) * i_rms * sqrt(2)
    sigma = t_lin / (2 * pi * d.r_gap ** 2 * d.stack * 1e-3)
    torque = t_lin / sqrt(1 + (sigma / (d.sigma_sat_kpa * 1e3)) ** 2)
    r_ph = rho_cu * d.n_phase * d.turn_len / d.cu_area
    w_m = d.rpm * 2 * pi / 60
    e_ph = d.kt_peak(t_mag) / 3 * 2 * w_m / sqrt(2)                  # rms phase back-EMF
    flow = mdot / p["rho"]
    return dict(j=j, i_rms=i_rms, torque=torque, t_lin=t_lin, sigma=sigma / 1e3,
                p_cu=p_cu, p_fe=p_fe, p_drv=p_drv, p_pump=dp_bar * 1e5 * flow / d.pump_eff, dp_bar=dp_bar,
                kt_ratio=torque / t_lin,
                t_max=t_max, t_mean=t_mean, t_mag=t_mag, t_out=t_coil_in + rise,
                t_skin=t_coil_in, v=v, re=re, h=h, flow_lpm=flow * 6e4,
                r_ph=r_ph, v_ph=e_ph + i_rms * r_ph, f_e=f_e,
                skin_mm=1e3 * sqrt(rho_cu / (pi * max(f_e, 1e-9) * MU0)))


def thermal_limit(d):
    """The current density at which the hottest copper reaches its limit with the
    pump at full pressure: how much heat the cooling can take, whether or not the
    iron can use the current."""
    lo, hi = 1e6, 400e6
    for _ in range(50):
        mid = sqrt(lo * hi)
        lo, hi = (mid, hi) if evaluate(d, mid)["t_max"] < d.t_cu_max else (lo, mid)
    return evaluate(d, lo)


def cooled(d, j):
    """Steady state at j with the pump pressure that gives the least total heat while
    keeping the copper under its limit and the flow under max_flow_lpm (more flow
    costs pump power but cools the copper, which lowers its resistance). None if
    that isn't enough."""
    best = None
    for dp in np.geomspace(1e-3, d.dp_bar, 28):
        r = evaluate(d, j, dp)
        if r["flow_lpm"] > d.max_flow_lpm:
            break
        if r["t_max"] <= d.t_cu_max and (best is None or total_loss(r) < total_loss(best)):
            best = r
    return best


def rated(d):
    """The continuous rating: the largest current that is both coolable and still
    worth having (torque per amp no more than kt_drop below its unsaturated value).
    Returns (result, which limit bound)."""
    j_th = thermal_limit(d)["j"]
    lo, hi = 1e5, j_th
    for _ in range(50):
        mid = sqrt(lo * hi)
        r = cooled(d, mid)
        ok = r is not None and r["kt_ratio"] > 1 - d.kt_drop
        lo, hi = (mid, hi) if ok else (lo, mid)
    r = cooled(d, lo)
    nxt = cooled(d, lo * 1.02)
    return r, ("cooling" if nxt is None else "saturation")


def total_loss(r):
    return r["p_cu"] + r["p_fe"] + r["p_drv"] + r["p_pump"]


def air_cooled(d, watts, fill=0.45, t_cu=120.0):
    """The same iron wound conventionally (solid wire) and limited to `watts` of
    copper loss: what the envelope does without the liquid."""
    cu = d.slot_area * fill / 2
    vol = cu * d.turn_len * d.slots
    rho = RHO_CU_20 * (1 + ALPHA_CU * (t_cu - 20))
    j = sqrt(watts / (rho * vol))
    amp_turns = j * cu                                               # per coil side
    t_lin = 3 * d.kw * (d.slots // 3) * d.b1(80) * d.r_gap * d.stack * 1e-3 * amp_turns * sqrt(2)
    sigma = t_lin / (2 * pi * d.r_gap ** 2 * d.stack * 1e-3)
    return dict(j=j, torque=t_lin / sqrt(1 + (sigma / (d.sigma_sat_kpa * 1e3)) ** 2), sigma=sigma / 1e3)


def magnet_grade(t):
    for grade, lim in (("N..SH", 150), ("N..UH", 180), ("N..EH", 200), ("SmCo", 300)):
        if t <= lim - 20:
            return grade
    return "none"


# ---- report -----------------------------------------------------------------------
def report(d, r):
    print(f"\n=== design point: Ø{d.od:.0f} x {d.stack + 2 * d.end_turn_mm:.0f} mm, "
          f"{'outrunner' if d.outrunner else 'stator outside, in the jacket'}, {d.slots}s/{d.poles}p ===")
    print(f"air-gap radius ...... {d.r_gap*1e3:.1f} mm; free bore Ø{d.bore_free_mm:.0f} mm; motor {d.mass():.2f} kg")
    side = sqrt(d.envelope_area) * 1e3
    print(f"conductor ........... {d.turns} turns/coil, {d.envelope_area*1e6:.1f} mm² envelope (~{side:.1f} mm square), "
          f"Ø{d.bore_d*1e3:.2f} mm bore, {d.cu_area*1e6:.1f} mm² copper")
    print(f"coolant ............. {COOLANTS[d.coolant]['label']}, in at {d.t_in:.0f} C, {r['dp_bar']:.2f} bar across each coil")
    print(f"  flow .............. {r['flow_lpm']:.2f} L/min, {r['v']:.1f} m/s in the bore, Re {r['re']:.0f}, "
          f"h {r['h']:.0f} W/m²K")
    print(f"  temperatures ...... skin {r['t_skin']:.0f} C, coolant out {r['t_out']:.0f} C, copper mean "
          f"{r['t_mean']:.0f} / max {r['t_max']:.0f} C, magnets ~{r['t_mag']:.0f} C ({magnet_grade(r['t_mag'])})")
    print(f"continuous .......... {r['j']/1e6:.1f} A/mm², {r['i_rms']:.0f} A rms per phase, "
          f"{r['v_ph']:.2f} V rms per phase at {d.rpm:.0f} rpm")
    print(f"  torque ............ {r['torque']:.1f} N·m ({r['t_lin']:.1f} before saturation; unsaturated shear {r['sigma']:.0f} kPa)")
    print(f"  torque density .... {r['torque']/d.mass():.1f} N·m/kg (motor only)")
    print(f"  losses ............ copper {r['p_cu']:.0f} W, drive + leads {r['p_drv']:.0f} W, iron {r['p_fe']:.0f} W, "
          f"pump {r['p_pump']:.0f} W  = {total_loss(r):.0f} W")
    print(f"  skin depth ........ {r['skin_mm']:.0f} mm at {r['f_e']:.0f} Hz, against a ~{(side - d.bore_d*1e3)/2:.1f} mm wall")


def main():
    d = Design()
    r, limit = rated(d)
    report(d, r)
    th = thermal_limit(d)
    print(f"\nrating is set by {limit}. With the pump at {d.dp_bar:.0f} bar and no flow cap the cooling would carry "
          f"{th['j']/1e6:.0f} A/mm² and {th['p_cu']/1e3:.0f} kW: {th['j']/r['j']:.1f}x the rated current, "
          f"for {th['torque']/r['torque']:.2f}x the torque")

    rs = SPECS["RS04"]
    i_rated = rs.tau_rated / rs.kt
    p_rated = 3 * i_rated ** 2 * rs.r_line / 2 * (1 + ALPHA_CU * 80)
    air = air_cooled(d, p_rated)
    print(f"\n=== baselines ===")
    print(f"RS04 (Ø{rs.od*1e3:.0f} mm, {rs.mass} kg with gearbox and driver): {rs.tau_rated/rs.ratio:.1f} N·m continuous "
          f"and {rs.tau_peak/rs.ratio:.1f} N·m peak at the motor, about {p_rated:.0f} W of copper loss at rated")
    print(f"this iron, solid wire, {p_rated:.0f} W: {air['torque']:.1f} N·m at {air['j']/1e6:.1f} A/mm² "
          f"-> liquid cooling gives {r['torque']/air['torque']:.1f}x the continuous torque for "
          f"{total_loss(r)/p_rated:.0f}x the heat")

    print(f"\n=== turns per coil, at the rated torque: current against bore size ===")
    print(f"  {'turns':>5} {'bore':>7} {'A rms':>7} {'V/ph':>6} {'N·m':>6} {'copper':>7} {'drive':>6} "
          f"{'pump':>6} {'total':>6} {'bar':>5} {'L/min':>6}")
    sweep_n = []
    for n in (1, 2, 3, 4, 6, 8, 12):
        dn = replace(d, turns=n)
        rn = cooled(dn, r["j"])
        if rn is None:
            print(f"  {n:5d} {dn.bore_d*1e3:6.2f}mm   can't be cooled at {d.dp_bar:.0f} bar")
            continue
        sweep_n.append((n, rn))
        print(f"  {n:5d} {dn.bore_d*1e3:6.2f}mm {rn['i_rms']:7.0f} {rn['v_ph']:6.2f} {rn['torque']:6.1f} "
              f"{rn['p_cu']:6.0f}W {rn['p_drv']:5.0f}W {rn['p_pump']:5.1f}W {total_loss(rn):5.0f}W "
              f"{rn['dp_bar']:5.2f} {rn['flow_lpm']:6.2f}")

    print(f"\n=== coolant return temperature, at the rated current (25 C ambient) ===")
    print(f"  {'in':>5} {'out':>5} {'skin':>5} {'N·m':>6} {'loss':>6} {'magnets':>8} {'grade':>6} {'radiator UA':>12}")
    sweep_t = []
    for t_in in (40, 60, 80, 100, 120, 140):
        dt = replace(d, t_in=t_in)
        rt = cooled(dt, r["j"])
        if rt is None:
            continue
        q = total_loss(rt) - rt["p_pump"]
        a, b = rt["t_out"] - 25, t_in - 25
        ua = q / ((a - b) / log(a / b))
        sweep_t.append((t_in, rt, ua))
        print(f"  {t_in:4.0f}C {rt['t_out']:4.0f}C {rt['t_skin']:4.0f}C {rt['torque']:6.1f} {q:5.0f}W {rt['t_mag']:7.0f}C "
              f"{magnet_grade(rt['t_mag']):>6} {ua:9.1f} W/K")

    print(f"\n=== variants (rated point of each) ===")
    print(f"  {'variant':42} {'N·m':>6} {'A/mm²':>6} {'loss':>6} {'N·m/kg':>7}  limit")
    variants = [
        ("design point", d),
        ("water-glycol (needs insulating breaks)", replace(d, coolant="glycol")),
        ("outrunner, stator on the hub", replace(d, outrunner=True)),
        ("copper limit 150 C", replace(d, t_cu_max=150.0)),
        ("saturation knee at 60 kPa", replace(d, sigma_sat_kpa=60.0)),
        ("saturation knee at 120 kPa", replace(d, sigma_sat_kpa=120.0)),
        ("no saturation at all", replace(d, sigma_sat_kpa=1e6)),
        ("stack 40 mm", replace(d, stack=40.0)),
    ]
    for label, dv in variants:
        rv, lim = rated(dv)
        print(f"  {label:42} {rv['torque']:6.1f} {rv['j']/1e6:6.1f} {total_loss(rv):5.0f}W {rv['torque']/dv.mass():7.1f}  {lim}")

    print(f"\n=== what it buys at the joint: gear ratio for RS04's output torque ===")
    print(f"  RS04: {rs.tau_rated} N·m continuous, {rs.tau_peak} N·m peak through {rs.ratio}:1, "
          f"about {p_rated:.0f} W to hold rated")
    print(f"  {'ratio':>5} {'cont. out':>10} {'heat to hold 35 N·m':>20} {'reflected inertia vs 9:1':>26}")
    js = np.geomspace(5e5, 2 * r["j"], 120)
    curve = [c for c in (cooled(d, j) for j in js) if c is not None]
    for ratio in (9, 6, 4, 3, 2, 1):
        need = rs.tau_rated / ratio / 0.96
        hit = next((c for c in curve if c["torque"] >= need and c["j"] <= r["j"] * 1.001), None)
        heat = f"{total_loss(hit):.0f} W" if hit else "over the rating"
        print(f"  {ratio:5d} {r['torque']*ratio*0.96:8.0f} N·m {heat:>20} {(ratio/9)**2:25.2f}x")

    pk = evaluate(d, 2 * r["j"])
    print(f"  short-term, at twice the rated current: {pk['torque']:.1f} N·m at the motor ({pk['p_cu']/1e3:.1f} kW of "
          f"copper loss); through 4:1 that is {pk['torque']*4*0.96:.0f} N·m against RS04's {rs.tau_peak}")
    plot(d, r, curve, sweep_n, sweep_t, air, p_rated)
    print(f"\nwrote {OUT / 'design.png'}")


def plot(d, r, curve, sweep_n, sweep_t, air, p_rated):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    loss = [total_loss(c) for c in curve]
    ax[0].plot(loss, [c["torque"] for c in curve], lw=2, label="liquid-cooled (this design)")
    ax[0].plot(loss, [c["t_lin"] for c in curve], "--", lw=1, color="gray", label="without saturation")
    ax[0].plot([p_rated], [air["torque"]], "o", color="tab:red", label="same iron, air-cooled budget")
    ax[0].plot([total_loss(r)], [r["torque"]], "s", color="k", label="rated point (torque per amp down 20%)")
    ax[0].set(xlabel="total heat (W)", ylabel="motor torque (N·m)", title="Torque against heat")
    ax[0].legend(fontsize=8)

    ns = [n for n, _ in sweep_n]
    for key, label in (("p_cu", "copper"), ("p_drv", "drive + leads"), ("p_pump", "pump")):
        ax[1].plot(ns, [r[key] / r["torque"] for _, r in sweep_n], "o-", label=label)
    ax[1].plot(ns, [total_loss(r) / r["torque"] for _, r in sweep_n], "k^-", label="total")
    ax[1].set(xlabel="turns per coil", ylabel="W per N·m at the rated torque", yscale="log",
              title="Turns: drive loss against bore size")
    ax[1].legend(fontsize=8)

    ts = [t for t, _, _ in sweep_t]
    ax[2].plot(ts, [total_loss(x) for _, x, _ in sweep_t], "o-", color="tab:blue")
    ax[2].set(xlabel="coolant return temperature (C)", ylabel="heat at the rated current (W)",
              title="Hotter loop: more loss, smaller radiator")
    ax[2].tick_params(axis="y", colors="tab:blue")
    ax2 = ax[2].twinx()
    ax2.plot(ts, [ua for _, _, ua in sweep_t], "s--", color="tab:orange")
    ax2.set_ylabel("radiator UA needed (W/K)", color="tab:orange")
    ax2.tick_params(axis="y", colors="tab:orange")
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "design.png", dpi=120)


if __name__ == "__main__":
    main()
