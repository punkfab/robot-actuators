"""
An RV (Rotate Vector) reducer actuator built around the reducer in Li et al. 2022
(Materials 15, 7850; see rcf.py), with its crankshaft's rolling contact fatigue life
carried through to the design.

The RV reducer, two stages:
  1. planetary: a sun on the motor shaft drives three planet gears, one on each crankshaft
  2. cycloid:   each crankshaft has two eccentrics (180 deg apart) that drive two cycloid
                discs through needle rollers running directly on the eccentric; the discs
                mesh with pins in the housing. The crankshafts ride in the carrier, which is
                the output: ratio i = 1 + (z_planet / z_sun) * z_pins

Taken from the paper (Table 1): 40 pins on R 82 mm, 39 lobes, e = 1.5 mm, three crankshafts
on r 46.77 mm, eccentric R 16.6 mm, rollers R 4 x 12 mm, rated output 800 N m, crankshaft at
585 rpm (so 14.6 rpm out), crankshaft 20CrNi2MoA carburized to a 1.0 mm case. Table 1 labels
the needle wheel 39 and the cycloid 40, but its short-amplitude coefficient k = e z / R_z =
0.7317 needs z = 40 pins, so the cycloid has 39 lobes (and it must have one fewer than the pins).

Chosen here (not in the paper): the planetary stage, pin and roller counts, disc thickness,
the motor, and the bearing clearance (fitted so the roller load reproduces the paper's
P_h = 1.56 GPa at 800 N m).

    python rv/design.py          (or: make rv)  -> out/design.json, out/design.png
"""

import json
import sys
from dataclasses import asdict, dataclass, replace
from math import floor, pi, sqrt, tau
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
sys.path.insert(0, str(HERE))
import rcf  # noqa: E402


@dataclass
class RV:
    # cycloid stage (the paper's)
    pins: int = 40
    r_pins: float = 82.0           # mm, pin circle radius
    pin_r: float = 3.5             # mm, pin radius (chosen: about half the 12.9 mm pitch)
    e: float = 1.5                 # mm, eccentricity
    disc_t: float = 12.0           # mm, each disc (= roller length, the paper's l1)
    n_crank: int = 3
    r_crank: float = 46.875        # mm (paper 46.77; set by the planetary centre distance below)
    ecc_r: float = 16.6            # mm, eccentric radius R2
    roller_r: float = 4.0          # mm, needle roller radius R1
    roller_gap: float = 1.2        # mm, between rollers (cage bars)
    clearance: float = 0.0         # mm, radial clearance of the needle bearing (fitted in main)
    roller_len: float = 0.0        # mm, effective roller length; 0 = disc_t (the paper's rollers span the disc)
    # planetary stage (chosen): module 1.25, 19 / 56 teeth -> centre distance 46.875 mm
    module: float = 1.25
    z_sun: int = 19
    z_planet: int = 56
    gear_w: float = 10.0

    @property
    def lobes(self): return self.pins - 1
    @property
    def k1(self): return self.e * self.pins / self.r_pins            # short amplitude coefficient, 0.7317
    @property
    def ratio(self): return 1 + self.z_planet / self.z_sun * self.pins
    @property
    def n_rollers(self):
        rp = self.ecc_r + self.roller_r
        return floor(tau * rp / (2 * self.roller_r + self.roller_gap))
    @property
    def centre_distance(self): return self.module * (self.z_sun + self.z_planet) / 2

    def pair(self):
        return rcf.Pair(r_roller=self.roller_r, r_ecc=self.ecc_r, length=self.roller_len or self.disc_t)


# ---- off-the-shelf bearings ------------------------------------------------------------------
@dataclass(frozen=True)
class Bearing:
    name: str
    d: float                       # mm bore
    D: float                       # mm outside
    B: float                       # mm width
    cr: float                      # N, basic dynamic radial rating
    c0r: float                     # N, basic static radial rating
    kind: str                      # roller | ball (life exponent 10/3 | 3)
    source: str

    def l10h(self, p, rpm):
        """basic rating life, hours, at equivalent load p (N) and speed rpm"""
        k = 10 / 3 if self.kind == "roller" else 3.0
        return (self.cr / p) ** k * 1e6 / (60 * rpm)


# Ratings from manufacturer datasheets or distributor listings (sources in README). Sizes are
# the standard ISO ones and all four are sold on Amazon under these designations.
# The eccentric bearing has no inner ring: its needles run on the crankshaft's eccentric, which
# must then be a raceway (>= 58 HRC, ground). The paper's carburized 670 HV surface is ~59 HRC.
# The outer ring presses into a plain bore in the disc, so the disc needs no hardening.
ECC = Bearing("RNA4906 needle roller bearing, no inner ring", 35, 47, 17, 25500, 39000, "roller",
              "RNA4906 (open) rating via Simply Bearings and Motion listings")
ECC_ALT = Bearing("HK3520 drawn-cup needle roller bearing", 35, 42, 20, 24800, 47500, "roller", "NTN HK3520 datasheet")
# RNA4906 internals, assumed (the catalogue pages found give no roller size): needles 3 mm on the
# 35 mm raceway, 13 mm effective length inside the 17 mm ring, 33 needles
RNA_DW, RNA_LW, RNA_Z = 3.0, 13.0, 33
JOURNAL = Bearing("32005 X tapered roller bearing", 25, 47, 15, 33000, 34000, "roller",
                  "Schaeffler 32005-X-XL datasheet")
MAIN = Bearing("61830 deep groove ball bearing", 150, 190, 20, 48800, 61000, "ball", "SKF 61830, via RS listing")
INPUT = Bearing("6004 deep groove ball bearing", 20, 42, 12, 9400, 5000, "ball", "NSK 6004, Amazon listing")
PIN = dict(name="DIN 6325 hardened dowel pin, m6", d=8.0, l=40.0)


def catalog_rv(rv):
    """the reducer with off-the-shelf parts: an RNA4906 in each disc (disc as wide as the bearing,
    a 35 mm hardened eccentric as the inner raceway), DIN 6325 8 mm pins"""
    rp = ECC.d / 2 + RNA_DW / 2
    gap = tau * rp / RNA_Z - RNA_DW                          # cage spacing that gives RNA_Z needles
    return replace(rv, disc_t=ECC.B, ecc_r=ECC.d / 2, pin_r=PIN["d"] / 2, roller_r=RNA_DW / 2,
                   roller_gap=gap - 1e-6, roller_len=RNA_LW, clearance=0.0)


# ---- loads ---------------------------------------------------------------------------------
def bearing_force(rv: RV, torque, phi):
    """resultant force on one crank needle bearing (N) at crank angle phi, paper eq. 1; torque N m.
    Uses z = pins in e*z (the paper's z5 = 40)."""
    t, ez, ro, k = torque * 1e3, rv.e * rv.pins, rv.r_crank, rv.k1
    return t / (2 * rv.n_crank * ez * ro) * np.sqrt(ez ** 2 + ro ** 2 + k ** 2 * ro ** 2
                                                    + 2 * ez * ro * np.cos(phi) - 2 * k * ez * ro * np.sin(phi))


def bearing_force_eq(rv: RV, torque):
    """equivalent bearing load over a revolution, paper eq. 2 (after Harris)"""
    t, ez, ro, k = torque * 1e3, rv.e * rv.pins, rv.r_crank, rv.k1
    num = ez ** 4 + 4 * (1 + k * k) * ez ** 2 * ro ** 2 + (1 + k * k) ** 2 * ro ** 4
    return t / 2 * (num / (rv.n_crank ** 4 * ez ** 4 * ro ** 4)) ** 0.25


def roller_loads(rv: RV, force):
    """load on each roller (N) of a caged needle bearing with radial clearance, under a radial
    force: Q = K delta^(10/9) per roller, delta = delta_r cos(psi) - c/2 (Harris), solved for
    delta_r. K from Palmgren's line-contact deflection, inner + outer contact in series."""
    z = rv.n_rollers
    psi = tau * np.arange(z) / z
    l = rv.roller_len or rv.disc_t
    # delta (mm) = 3.84e-5 Q^0.9 / l^0.8 per contact -> two contacts: Q = K delta^(1/0.9)
    k = (1 / (2 * 3.84e-5 / l ** 0.8)) ** (1 / 0.9)
    lo, hi = 0.0, 0.1
    for _ in range(100):
        d = 0.5 * (lo + hi)
        q = k * np.maximum(d * np.cos(psi) - rv.clearance / 2, 0) ** (1 / 0.9)
        if np.sum(q * np.cos(psi)) > force:
            hi = d
        else:
            lo = d
    return q


def p_h(rv: RV, torque):
    """peak contact pressure on the eccentric (MPa) at the peak bearing force over a revolution"""
    f = float(np.max(bearing_force(rv, torque, np.linspace(0, tau, 721))))
    q = float(np.max(roller_loads(rv, f)))
    return rv.pair().p0_of_load(q), q, f


def fit_clearance(rv: RV, torque=800.0, target=1560.0):
    """radial clearance that reproduces the paper's P_h at the rated torque"""
    lo, hi = 0.0, 0.2
    for _ in range(60):
        c = 0.5 * (lo + hi)
        if p_h(replace(rv, clearance=c), torque)[0] > target:
            hi = c
        else:
            lo = c
    return 0.5 * (lo + hi)


# ---- cycloid profile and the pin contact -------------------------------------------------------
def disc_profile(rv: RV, steps=2000, mod=0.0):
    """disc outline in its own frame (pin-wheel cycloid, offset by the pin radius + mod), mm"""
    R, Rr, E, N = rv.r_pins, rv.pin_r + mod, rv.e, rv.pins
    t = tau * np.arange(steps) / steps
    psi = np.arctan2(np.sin((1 - N) * t), R / (E * N) - np.cos((1 - N) * t))
    x = R * np.cos(t) - Rr * np.cos(t + psi) - E * np.cos(N * t)
    y = -R * np.sin(t) + Rr * np.sin(t + psi) + E * np.sin(N * t)
    return x, y


def curvature(x, y):
    """signed curvature of a closed curve (periodic central differences)"""
    d = lambda v: (np.roll(v, -1) - np.roll(v, 1)) / 2
    dx, dy = d(x), d(y)
    ddx, ddy = np.roll(x, -1) - 2 * x + np.roll(x, 1), np.roll(y, -1) - 2 * y + np.roll(y, 1)
    return (dx * ddy - dy * ddx) / (dx * dx + dy * dy) ** 1.5


def undercut_ok(rv: RV):
    """the actual profile must not loop: its tangent never reverses and it never crosses
    itself. (The pin-centre curve's own radius at the tooth valleys is tiny, 0.2 mm here, but
    the inward offset by the pin radius opens it out, so that radius is not the criterion.)
    Returns the tightest radius of the actual profile and the verdict."""
    from shapely.geometry import LinearRing
    x, y = disc_profile(rv, steps=40000)
    tx, ty = np.roll(x, -1) - x, np.roll(y, -1) - y
    cos_turn = (tx * np.roll(tx, -1) + ty * np.roll(ty, -1)) / (np.hypot(tx, ty) * np.hypot(np.roll(tx, -1), np.roll(ty, -1)))
    ok = bool(cos_turn.min() > 0 and LinearRing(np.c_[x, y]).is_simple)
    return float(np.min(1 / np.abs(curvature(x, y)))), ok


def pin_contact(rv: RV, torque, share=0.55):
    """peak pin force and contact pressure. Pin force by the classic estimate for a cycloid with
    about half the pins in contact: F_max = 4 T_disc / (K1 z_c R_p) (e.g. Blanche & Yang 1989;
    Li et al.), T_disc = share x T (two discs, unequal sharing). Contact: pin against the
    profile at its tightest curvature in the working (concave) zone, line contact over the disc
    thickness."""
    t_disc = share * torque * 1e3
    f = 4 * t_disc / (rv.k1 * rv.lobes * rv.r_pins)
    x, y = disc_profile(rv, steps=8000)
    k = curvature(x, y)
    # working flank: the disc's concave-to-pin regions are where |rho| is largest; take the
    # convex lobe tip radius (worst case) for the profile radius
    rho = float(np.min(1 / np.abs(k)))
    e_star = 1 / (2 * (1 - 0.3 ** 2) / 210e3)
    r_eq = 1 / (1 / rv.pin_r + 1 / rho)
    p0 = sqrt(f / rv.disc_t * e_star / (pi * r_eq))
    return f, p0, rho


# ---- motor ---------------------------------------------------------------------------------------
def motor_for(rv: RV, t_out_rated, t_out_peak, n_out_max, eta=0.85, bore=20.0):
    """smallest stack of a 12s/14p motor (rollerscrew's model) that gives the rated input torque
    continuously and the peak for a few seconds, at 48 V up to the top speed"""
    import importlib.util
    spec = importlib.util.spec_from_file_location("rollerscrew_design", HERE.parent / "rollerscrew" / "design.py")
    rsd = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rsd)
    M = rsd.Motor                                  # the roller-screw actuator's motor model
    t_rated = t_out_rated / rv.ratio / eta
    t_peak = t_out_peak / rv.ratio / eta
    w_max = n_out_max * rv.ratio * tau / 60
    best = None
    for depth in (8.0, 10.0, 12.0, 14.0):
        for stack in np.arange(10, 90, 2.0):
            for bore_x in (bore, bore + 10, bore + 20, bore + 30):
                m = M(bore=bore_x, stack=float(stack), slot_depth=depth, h_eff=25.0, r_int=0.3)
                if m.torque(m.j_cont) < t_rated or m.torque(m.j_peak) < t_peak:
                    continue
                # turns: the largest that still reaches w_max at the rated torque
                ok = None
                for n in range(1, 80):
                    mm = replace(m, turns=n)
                    j = t_rated / mm.torque_lin(1.0)
                    if mm.speed(j) >= w_max:
                        ok = mm
                if ok is None:
                    continue
                mass = sum(ok.masses().values())
                if best is None or mass < best[0]:
                    best = (mass, ok)
                break
    return best[1], t_rated, t_peak, w_max


# ---- service life of the crankshaft eccentric ---------------------------------------------------
# rcf.py's elastic model runs 4.7x short of the paper's FEM life (13% in damage): calibrate once,
# on the paper's ICRS case. Residual-stress states the elastic model cannot see enter as the
# paper's own FEM life ratios (Fig. 24), relative to "none".
PAPER_RS_FACTOR = {k: v / rcf.PAPER["none"] for k, v in rcf.PAPER.items()}


def crank_life(rv: RV, torque, case: rcf.Case, crank_rpm, calib):
    """crack-initiation life of the eccentric at a steady torque: cycles (as the paper counts
    them, one per crank revolution) and hours"""
    p0, q, f = p_h(rv, torque)
    nf, depth, _ = rcf.life(rv.pair(), replace(case, residual="none"), p0,
                            depths=np.arange(0.01, 0.31, 0.01), plane_step=10.0)
    nf *= calib * PAPER_RS_FACTOR[case.residual]
    return dict(p0=p0, roller_load=q, bearing_force=f, nf=nf, depth=depth, hours=nf / (crank_rpm * 60))


def calibrate():
    """life factor that brings rcf.py's elastic 'none' result to the paper's FEM value at 1.56 GPa"""
    nf, _, _ = rcf.life(rcf.Pair(), rcf.Case(residual="none"), 1560.0, depths=np.arange(0.01, 0.31, 0.01), plane_step=10.0)
    return rcf.PAPER["none"] / nf, nf


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(exist_ok=True)
    rv = RV()
    print(f"RV: {rv.pins} pins on R{rv.r_pins}, {rv.lobes} lobes, e {rv.e}, K1 {rv.k1:.4f} (paper 0.7317); "
          f"planetary {rv.z_sun}/{rv.z_planet} m{rv.module} -> centre distance {rv.centre_distance:.3f} mm (paper 46.77)")
    print(f"ratio 1 + {rv.z_planet}/{rv.z_sun} x {rv.pins} = {rv.ratio:.2f}; crank turns {rv.pins}x output relative to the carrier")
    rho, ok = undercut_ok(rv)
    print(f"profile: pin radius {rv.pin_r} mm, tightest profile radius {rho:.2f} mm, "
          f"{'no loops (no undercut)' if ok else 'LOOPS: undercut'}")
    print(f"needle bearing: {rv.n_rollers} rollers Ø{2 * rv.roller_r} x {rv.disc_t} on the R{rv.ecc_r} eccentric")

    t_rated = 800.0
    c_fit = fit_clearance(rv, t_rated)
    rv = replace(rv, clearance=round(c_fit, 4))
    p0, q, f = p_h(rv, t_rated)
    f_eq = bearing_force_eq(rv, t_rated)
    print(f"\nat {t_rated:.0f} N m: bearing force {f:.0f} N peak (paper Fig. 7 about 5200), "
          f"{f_eq:.0f} N equivalent (eq. 2)")
    q0 = float(np.max(roller_loads(replace(rv, clearance=0.0), f)))
    print(f"  zero clearance: max roller {q0:.0f} N -> P_h {rv.pair().p0_of_load(q0) / 1e3:.2f} GPa; "
          f"radial clearance {rv.clearance * 1e3:.1f} um reproduces the paper's 1.56 GPa (max roller {q:.0f} N)")
    fp, pp, rho_p = pin_contact(rv, t_rated)
    print(f"  pins: peak pin force {fp:.0f} N, contact {pp / 1e3:.2f} GPa (profile radius {rho_p:.1f} mm)")

    calib, nf_raw = calibrate()
    print(f"\nfatigue model calibration: elastic 'none' life at 1.56 GPa {nf_raw:.3g} -> x{calib:.2f} to the paper's "
          f"{rcf.PAPER['none']:.3g}; residual-stress factors from the paper: "
          + ", ".join(f"{k} {v:.2f}" for k, v in PAPER_RS_FACTOR.items()))
    crank_rpm = 585.0
    base = crank_life(rv, t_rated, rcf.Case(residual="icrs"), crank_rpm, calib)
    print(f"  as carburized, 800 N m steady at {crank_rpm:.0f} rpm: N_f {base['nf']:.3g} = {base['hours']:.0f} h "
          f"(paper test: 2150 h)")

    # load-life: the eccentric at other torques, and the process choices
    torques = np.array([400, 500, 600, 700, 800, 900, 1000, 1200])
    curves = {}
    for rs in ("sprs", "icrs", "none", "itrs"):
        curves[rs] = [crank_life(rv, float(t), rcf.Case(residual=rs), crank_rpm, calib)["hours"] for t in torques]
    print("\nhours at steady torque, crank at 585 rpm:")
    print("  torque " + "".join(f"{t:>9.0f}" for t in torques))
    for rs, h in curves.items():
        print(f"  {rs:6} " + "".join(f"{x:9.3g}" for x in h))
    sens = {}
    for name, c in (("surface 620 HV", rcf.Case(hv_surface=620.0)), ("surface 700 HV", rcf.Case(hv_surface=700.0)),
                    ("case 0.6 mm", rcf.Case(chd=0.6)), ("case 1.4 mm", rcf.Case(chd=1.4))):
        sens[name] = crank_life(rv, t_rated, c, crank_rpm, calib)["hours"]
    print("  case and hardness at 800 N m (ICRS): " + ", ".join(f"{k} {v:.0f} h" for k, v in sens.items()))

    # the same reducer with off-the-shelf bearings and pins
    cat = catalog_rv(rv)
    rho_c, ok_c = undercut_ok(cat)
    f_eq_c = bearing_force_eq(cat, t_rated)
    f_pk_c = float(np.max(bearing_force(cat, 2 * t_rated, np.linspace(0, tau, 721))))
    n_in = crank_rpm * RV().z_planet / RV().z_sun * (1 - 1 / rv.ratio) ** -1
    fp_c, pp_c, _ = pin_contact(cat, t_rated)
    lives = dict(
        eccentric=ECC.l10h(f_eq_c, crank_rpm),
        journal=JOURNAL.l10h(f_eq_c, crank_rpm),                        # assumed: each journal carries about one eccentric's load
    )
    print(f"\noff-the-shelf version: discs {cat.disc_t:.0f} mm (the RNA4906 width), hardened eccentric Ø{2 * cat.ecc_r:.0f}, "
          f"pins Ø{2 * cat.pin_r:.0f} x {PIN['l']:.0f} ({PIN['name']}); profile {'no loops' if ok_c else 'LOOPS'} "
          f"(tightest radius {rho_c:.2f} mm)")
    print(f"  eccentric bearing {ECC.name}: P = {f_eq_c:.0f} N (eq. 2) at {crank_rpm:.0f} rpm -> L10 {lives['eccentric']:.0f} h; "
          f"peak force at 2x rated torque {f_pk_c:.0f} N = C0r / {ECC.c0r / f_pk_c:.1f}")
    print(f"  crank journals {JOURNAL.name}: L10 {lives['journal']:.0f} h at the same load (assumption)")
    print(f"  main bearings {MAIN.name}: Cr {MAIN.cr / 1e3:.1f} kN each; under pure torque they carry almost nothing "
          f"(three cranks balance), their load is the arm's external moment (not modelled), output {crank_rpm / rv.pins:.1f} rpm")
    print(f"  input shaft {INPUT.name}: three planets balance the sun's radial load; input {n_in:.0f} rpm at rated")
    print(f"  pins: {fp_c:.0f} N peak pin force, {pp_c / 1e3:.2f} GPa")
    print(f"  compare the paper's integral needle bearing (rollers on the crank): {base['hours']:.0f} h as carburized")
    # the paper's model on this eccentric: the RNA4906's needles on the carburized crankshaft
    p0_c, q_c, _ = p_h(cat, t_rated)
    rcf_c = crank_life(cat, t_rated, rcf.Case(residual="icrs"), crank_rpm, calib)
    lives["eccentric_rcf"] = rcf_c["hours"]
    print(f"  the paper's crankshaft model on the RNA4906 eccentric ({cat.n_rollers} needles Ø{2 * cat.roller_r:.0f} x "
          f"{cat.roller_len:.0f}, assumed): max needle {q_c:.0f} N, P_h {p0_c / 1e3:.2f} GPa (paper's design 1.56), "
          f"crack initiation {rcf_c['hours']:.3g} h as carburized")

    # actuator: motor on the input
    t_peak, n_out_max = 1600.0, 25.0                     # start/stop torque 2x rated; top speed
    mot, t_in, t_in_pk, w_max = motor_for(rv, t_rated, t_peak, n_out_max)
    ph = mot.phase()
    print(f"\nmotor: input {t_in:.2f} N m rated, {t_in_pk:.2f} N m peak, up to {w_max * 60 / tau:.0f} rpm "
          f"(output {n_out_max:.0f} rpm); chose Ø{mot.od:.0f} x {mot.stack:.0f} mm stack, slots {mot.slot_depth:.0f} mm deep, "
          f"bore {mot.bore:.0f}, {mot.turns} turns: {mot.torque(mot.j_cont):.1f} N m continuous, {mot.torque(mot.j_peak):.1f} peak, "
          f"Kt {ph['kt']:.3f} N m/A, {ph['i_per_j'] * mot.j_peak:.0f} A peak, {sum(mot.masses().values()):.2f} kg active")

    fig, ax = plt.subplots(1, 3, figsize=(12.5, 3.8))
    phi = np.linspace(0, 3 * pi, 541)
    for k in range(rv.n_crank):
        ax[0].plot(np.degrees(phi), bearing_force(rv, t_rated, phi + k * tau / 3), lw=1.3, label=f"crank bearing {k + 1}")
    ax[0].axhline(f_eq, color="0.5", lw=0.8, ls="--", label="equivalent (eq. 2)")
    ax[0].set(xlabel="crank angle (deg)", ylabel="bearing force (N)", title="needle-bearing load at 800 N m (paper Fig. 7)")
    ax[0].legend(frameon=False, fontsize=7.5)
    ql = roller_loads(rv, f)
    ax[1].bar(np.arange(rv.n_rollers), ql, color="#1f4e79", label=f"clearance {rv.clearance * 1e3:.0f} um")
    ax[1].plot(np.arange(rv.n_rollers), roller_loads(replace(rv, clearance=0.0), f), "o", color="#c87a2f", ms=4, label="no clearance")
    ax[1].set(xlabel="roller (0 = load line)", ylabel="roller load (N)", title=f"{rv.n_rollers} rollers at the peak force")
    ax[1].legend(frameon=False, fontsize=8)
    for rs, col in (("sprs", "#b03a2e"), ("icrs", "#c87a2f"), ("none", "0.5"), ("itrs", "#1f4e79")):
        ax[2].loglog(torques, curves[rs], "o-", ms=3, color=col, label=rs.upper())
    ax[2].axhline(2150, color="k", lw=0.6, ls=":")
    ax[2].text(410, 2300, "paper's test, 2150 h", fontsize=7.5)
    ax[2].set(xlabel="output torque (N m)", ylabel="eccentric life (h at 585 rpm)", title="crankshaft load-life by surface state")
    ax[2].legend(frameon=False, fontsize=8)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "design.png", dpi=150)
    out = dict(rv=asdict(cat), paper_rv=asdict(rv), ratio=rv.ratio,
               bearings={k: asdict(b) for k, b in dict(eccentric=ECC, eccentric_alt=ECC_ALT, journal=JOURNAL, main=MAIN,
                                                        input=INPUT).items()},
               pin=PIN, catalog_lives=lives, catalog_peak_force=f_pk_c, catalog_pin=dict(force=fp_c, p0=pp_c), n_rollers=rv.n_rollers, k1=rv.k1, undercut_rho=rho,
               rated=dict(torque=t_rated, bearing_peak=f, bearing_eq=f_eq, roller_peak=q, p_h=p0, pin_force=fp, pin_p0=pp),
               calibration=calib, life_rated=base, life_curves=dict(torque=torques.tolist(), **curves), sensitivity=sens,
               motor=dict(od=mot.od, stack=mot.stack, slot_depth=mot.slot_depth, bore=mot.bore, turns=mot.turns, kt=ph["kt"],
                          t_cont=mot.torque(mot.j_cont), t_peak=mot.torque(mot.j_peak), i_peak=ph["i_per_j"] * mot.j_peak,
                          r_phase=ph["r"], mass=sum(mot.masses().values()), t_in_rated=t_in, t_in_peak=t_in_pk, w_max=w_max))
    (OUT / "design.json").write_text(json.dumps(out, indent=1, default=float))
    print(f"\nwrote {OUT / 'design.json'}, {OUT / 'design.png'}")


if __name__ == "__main__":
    main()
