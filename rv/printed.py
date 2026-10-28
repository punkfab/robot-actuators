"""
A 3D-printable RV reducer: the same architecture as design.py / cad.py (sun -> three planets
on crankshafts -> two cycloid discs 180 deg apart on fixed pins -> carrier output), but every
part printed except standard deep-groove ball bearings and hardened steel dowel pins.

  bought      6802-2RS x12 (two per eccentric), 6800-2RS x6 (crank journals), 6814-2RS x2 (main),
              DIN 6325 dowels: 5 x 24 x30 (pin ring), 5 x 45 x3 (carrier pillars), 5 x 45 x3
              (crankshaft cores), a NEMA 17 stepper
  printed     housing, input cover, carrier plates (x2), discs (x2), crankshafts (two pieces each,
              hex-keyed, on the steel core), sun, planets (on the crankshafts)

The design question for a printed reducer is which printed feature gives out first, so
analyse() ranks them: pin-on-lobe contact, gear tooth bending, crank hex key, against
allowable stresses for the print material (assumed, stated below), and the ball bearings
against their catalogue ratings.

    python rv/printed.py          (or: make rv-printed)  -> out/printed/
"""

import json
import sys
from dataclasses import dataclass
from math import cos, pi, radians, sin, sqrt, tau
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "printed"
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402

# ---- geometry ------------------------------------------------------------------------------------
P = dz.RV(pins=30, r_pins=45.0, pin_r=2.5, e=1.0, disc_t=10.0, n_crank=3, r_crank=22.0, module=1.0,
          z_sun=12, z_planet=32, gear_w=6.0)
MOD = 0.10                      # mm, disc profile offset: print tolerance + running clearance
GROOVE_CLEAR = 0.08             # mm on the radius: the pins float in their grooves and roll on the lobes
GAP, GAP_F = 1.0, 1.5           # disc to disc, plate to disc: 2 x 10 + 1 + 2 x 1.5 = 24 = the pin length
PIN_L = 2 * P.disc_t + GAP + 2 * GAP_F
PLATE_T = 10.5                  # carrier plates: 2 x 10.5 + 24 = 45 = the pillar dowel length
R_PLATE = 35.0                  # = main bearing bore
R_PILLAR, R_CORE = 2.5, 2.5     # Ø5 dowels
CRANK0 = 90.0                   # as cad.py: lets disc 2 be disc 1 flipped over

# bearings: (d, D, B, Cr N) -- ratings are listing values (sources in README)
B_ECC = (15.0, 24.0, 5.0, 1900.0)      # 6802-2RS, two side by side per disc
B_JRN = (10.0, 19.0, 5.0, 1800.0)      # 6800-2RS
B_MAIN = (70.0, 90.0, 10.0, 12100.0)   # 6814-2RS
NEMA = dict(square=42.3, boss=22.0, bolt_spacing=31.0, shaft=5.0, flat=0.5, body=40.0)

# print material (PETG, assumed allowables for repeated loading; change for PLA / nylon / CF)
MAT = dict(name="PETG", E=2100.0, nu=0.38, contact=40.0, bending=20.0, shear=15.0)
T_MOTOR = 0.40                  # N m: a typical NEMA 17 at modest speed (assumed; holding torque is higher)
ETA = 0.75                      # assumed for a printed reducer


def ratio():
    return P.ratio


# ---- analysis ------------------------------------------------------------------------------------
def pin_force(t_out, share=0.55):
    """peak pin force (N): 4 T_disc / (K1 z_c R_p), as design.pin_contact"""
    return 4 * share * t_out * 1e3 / (P.k1 * P.lobes * P.r_pins)


def lobe_contact(t_out):
    """steel pin on a printed lobe, line contact over the disc thickness (MPa)"""
    x, y = dz.disc_profile(P, steps=8000, mod=MOD)
    rho = float(np.min(1 / np.abs(dz.curvature(x, y))))
    e_star = 1 / ((1 - 0.3 ** 2) / 210e3 + (1 - MAT["nu"] ** 2) / MAT["E"])
    r_eq = 1 / (1 / P.pin_r + 1 / rho)
    return sqrt(pin_force(t_out) / P.disc_t * e_star / (pi * r_eq)), rho


def tooth_bending(t_out):
    """Lewis bending at the sun's root (MPa): input torque over three planet meshes"""
    t_in = t_out / P.ratio / ETA
    ft = t_in * 1e3 / (P.module * P.z_sun / 2) / P.n_crank
    y_lewis = 0.245                                           # z = 12, 20 deg full depth
    return ft / (P.gear_w * P.module * y_lewis)


def hex_shear(t_out):
    """the hex key joining the two crankshaft pieces (MPa): crank torque on an 8 mm hex, 8 mm deep"""
    t_crank = t_out / P.ratio / ETA * P.z_planet / P.z_sun / P.n_crank
    af, depth = 8.0, 8.0
    # bearing stress on the hex flats ~ T / (n_flats x lever x flat area)
    return t_crank * 1e3 / (6 * (af / 2) * (af / sqrt(3) * depth))


def torque_limit(fn, allow, lo=0.1, hi=500.0):
    for _ in range(60):
        mid = sqrt(lo * hi)
        if fn(mid) > allow:
            hi = mid
        else:
            lo = mid
    return lo


# ---- ring-pin options: efficiency, using cycloidal/efficiency.py's contact model ----------------
RING_OPTIONS = [
    ("printed teeth, part of the housing", "printed", "lobe slides on printed plastic"),
    ("steel dowels, fixed, dry", "steel_dry", "lobe slides on steel"),
    ("steel dowels, fixed, greased", "steel_lubed", "lobe slides on greased steel"),
    ("steel dowels, floating in their grooves", "sleeve", "pin rolls on the lobe, slides in its groove"),
    ("needle-bearing rollers on posts", "needle", "rolling throughout"),
]


def ring_options():
    """asymptotic efficiency for each way of making the ring pins. The RV takes its output
    through the crankshafts' ball bearings, not sliding output pins, so the output term is a
    rolling one; the printed planetary stage is in series. Absolute values are the model's
    uncalibrated estimates; the ranking is what the friction ratios establish."""
    sys.path.insert(0, str(HERE.parent / "cycloidal"))
    import efficiency as ef
    eta_pl = 1 - ef.PLANET_STAGE_LOSS["printed"]
    rows = []
    for name, contact, how in RING_OPTIONS:
        eta_c = ef.predict(contact, "needle")
        eta = eta_pl * eta_c
        rows.append(dict(name=name, contact=contact, how=how, mu=ef.CONTACT_FRICTION[contact], eta_cyclo=eta_c, eta=eta,
                         eta_back=max(0.0, 2 - 1 / eta),                  # back-driving, for a simple Coulomb loss
                         t_out=T_MOTOR * P.ratio * eta))
    # contact pressure: a printed tooth is as soft as the lobe, which halves E* (pressure / sqrt 2)
    e_ps = 1 / ((1 - 0.3 ** 2) / 210e3 + (1 - MAT["nu"] ** 2) / MAT["E"])
    e_pp = 1 / (2 * (1 - MAT["nu"] ** 2) / MAT["E"])
    return rows, dict(e_star_steel_pin=e_ps, e_star_printed_tooth=e_pp, pressure_ratio=sqrt(e_pp / e_ps))


def analyse():
    t_nominal = T_MOTOR * P.ratio * ETA
    rho_ok = dz.undercut_ok(P)
    limits = {
        "pin on printed lobe (contact)": torque_limit(lambda t: lobe_contact(t)[0], MAT["contact"]),
        "sun tooth root (bending)": torque_limit(tooth_bending, MAT["bending"]),
        "crankshaft hex key (bearing)": torque_limit(hex_shear, MAT["shear"]),
    }
    f_eq = dz.bearing_force_eq(P, t_nominal)
    crank_rpm = 60.0                                     # an example: 60 rpm crank = 1.5 rpm out
    l10 = lambda cr, p, n_brg=1: (cr * n_brg / p) ** 3 * 1e6 / (60 * crank_rpm)
    lives = {"eccentric, 2 x 6802-2RS": l10(B_ECC[3], f_eq, 2), "journal, 6800-2RS": l10(B_JRN[3], f_eq)}
    return dict(ratio=P.ratio, k1=P.k1, profile_ok=rho_ok[1], tightest_radius=rho_ok[0], t_nominal=t_nominal,
                limits=limits, governing=min(limits, key=limits.get), bearing_force_eq=f_eq, lives_h=lives,
                lobe_p0_nominal=lobe_contact(t_nominal)[0], crank_rpm=crank_rpm)


# ---- CAD -------------------------------------------------------------------------------------------
def build():
    from build123d import (Align, Box, BuildLine, BuildPart, BuildSketch, Cylinder, Kind, Plane, Pos,
                           RegularPolygon, Rot, Solid, Spline, extrude, make_face, mirror, offset)
    import cad as steel_cad                                 # spur_gear and the shared helpers

    def cyl(r, h, z0=0.0, x=0.0, y=0.0):
        return Pos(x, y, z0) * Cylinder(r, h, align=(Align.CENTER, Align.CENTER, Align.MIN))

    def ring(r_in, r_out, z0, h, x=0.0, y=0.0):
        c = cyl(r_out, h, z0, x, y)
        return c - cyl(r_in, h + 0.2, z0 - 0.1, x, y) if r_in > 0 else c

    def crank_xy(k):
        a = radians(CRANK0) + tau * k / P.n_crank
        return P.r_crank * cos(a), P.r_crank * sin(a)

    def pillar_xy(k):
        a = radians(CRANK0) + tau * (k + 0.5) / P.n_crank
        return P.r_crank * cos(a), P.r_crank * sin(a)

    # stations (z, mm); the motor face at z = 0, output at the far end
    Z = dict(cover=0.0, gears=10.0)
    Z["plate_r"] = Z["gears"] + P.gear_w + 3.0
    Z["disc1"] = Z["plate_r"] + PLATE_T + GAP_F
    Z["disc2"] = Z["disc1"] + P.disc_t + GAP
    Z["plate_f"] = Z["disc2"] + P.disc_t + GAP_F
    Z["end"] = Z["plate_f"] + PLATE_T
    Z["pins"] = Z["plate_r"] + PLATE_T

    E = P.e
    p = {}
    # pins (steel) and housing (printed): pins half-sunk in grooves of the pin case
    p["pin"] = [cyl(P.pin_r, PIN_L, Z["pins"], P.r_pins * cos(tau * k / P.pins), P.r_pins * sin(tau * k / P.pins))
                for k in range(P.pins)]
    r_out = B_MAIN[1] / 2 + 5.0
    h = ring(B_MAIN[1] / 2, r_out, Z["plate_r"], Z["end"] - Z["plate_r"])
    h = h.fuse(ring(P.r_pins, r_out - 1, Z["pins"], PIN_L))
    h = h.fuse(ring(B_MAIN[1] / 2 - 4.0, r_out, Z["end"], 2.0))                     # output-side lip retains the bearing
    h = h.cut(*[cyl(P.pin_r + GROOVE_CLEAR, PIN_L + 0.2, Z["pins"] - 0.1, P.r_pins * cos(tau * k / P.pins), P.r_pins * sin(tau * k / P.pins))
                for k in range(P.pins)])
    p["housing"] = [h]

    # discs: one part, disc 2 flipped
    x, y = dz.disc_profile(P, steps=P.lobes * 64, mod=MOD)
    with BuildPart() as bp:
        with BuildSketch():
            with BuildLine():
                Spline(*zip(x.tolist(), y.tolist()), periodic=True)
            make_face()
        extrude(amount=P.disc_t)
    d = bp.part
    holes = [cyl(B_ECC[1] / 2, P.disc_t + 0.2, -0.1, *crank_xy(k)) for k in range(P.n_crank)]
    holes += [cyl(R_PILLAR + E + 0.5, P.disc_t + 0.2, -0.1, *pillar_xy(k)) for k in range(P.n_crank)]
    holes.append(cyl(8.0, P.disc_t + 0.2, -0.1))
    d = d.cut(*holes)
    p["disc"] = [Pos(E, 0, Z["disc1"]) * d, Pos(-E, 0, Z["disc2"]) * mirror(d, about=Plane.YZ)]

    # eccentric bearings: two 6802 per disc
    p["ecc_bearing"] = [ring(B_ECC[0] / 2, B_ECC[1] / 2, z0 + j * B_ECC[2], B_ECC[2], crank_xy(k)[0] + s * E, crank_xy(k)[1])
                        for k in range(P.n_crank) for s, z0 in ((1, Z["disc1"]), (-1, Z["disc2"])) for j in (0, 1)]

    # carrier plates (printed) with journal and pillar seats; pillars (steel dowels)
    def plate(z0):
        f = ring(0.0, R_PLATE, z0, PLATE_T)
        cuts = [cyl(B_JRN[1] / 2, PLATE_T + 0.2, z0 - 0.1, *crank_xy(k)) for k in range(P.n_crank)]
        cuts += [cyl(R_PILLAR, PLATE_T + 0.2, z0 - 0.1, *pillar_xy(k)) for k in range(P.n_crank)]
        cuts.append(cyl(6.0, PLATE_T + 0.2, z0 - 0.1))
        return f.cut(*cuts)

    p["carrier_rear"] = [plate(Z["plate_r"])]
    out_plate = plate(Z["plate_f"])
    out_plate = out_plate.cut(*[cyl(1.25, 8.0, Z["end"] - 8.0, 28.0 * cos(a), 28.0 * sin(a))   # M3 tapping holes, output
                                for a in np.radians(np.arange(15, 360, 60))])
    p["carrier_front"] = [out_plate]
    p["pillar"] = [cyl(R_PILLAR, 2 * PLATE_T + PIN_L, Z["plate_r"], *pillar_xy(k)) for k in range(P.n_crank)]
    p["journal_bearing"] = [ring(B_JRN[0] / 2, B_JRN[1] / 2, z0 + (PLATE_T - B_JRN[2]) / 2, B_JRN[2], *crank_xy(k))
                            for k in range(P.n_crank) for z0 in (Z["plate_r"], Z["plate_f"])]
    p["main_bearing"] = [ring(B_MAIN[0] / 2, B_MAIN[1] / 2, Z["plate_r"] + PLATE_T - B_MAIN[2], B_MAIN[2]),
                         ring(B_MAIN[0] / 2, B_MAIN[1] / 2, Z["plate_f"], B_MAIN[2])]

    # crankshafts: piece A = planet + rear journal + hex boss; piece B = both eccentrics + front journal
    # with a hex socket; a steel core dowel through both
    pitch = 360.0 / P.z_planet
    sun = steel_cad.spur_gear(P.module, P.z_sun, P.gear_w + 2.0, bore=0)
    from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
    sun = offset(sun, amount=-0.15, kind=Kind.INTERSECTION)                         # printed-gear backlash
    sun = Solid(BRepBuilderAPI_NurbsConvert(sun.solids()[0].wrapped, True).Shape())
    shaft_d = NEMA["shaft"]
    sun = sun - cyl(shaft_d / 2 + 0.1, P.gear_w + 4, -1.0)
    sun = sun + Pos(0, shaft_d / 2 - NEMA["flat"] + 1.5, P.gear_w / 2 + 1) * Box(shaft_d, 3.0, P.gear_w + 2.0)   # the D flat
    sun = Pos(0, 0, Z["gears"] - 1.0) * sun
    p["sun"] = [sun]
    hexr = 8.0 / sqrt(3)
    pieces_a, pieces_b, cores = [], [], []
    shared = steel_cad.shared
    for k in range(P.n_crank):
        cx, cy = crank_xy(k)
        planet0 = steel_cad.spur_gear(P.module, P.z_planet, P.gear_w, bore=0)
        best = None
        for i in range(24):
            t = pitch * i / 24
            g = Pos(cx, cy, Z["gears"]) * Rot(0, 0, t) * planet0
            v = shared(sun, g)
            if best is None or v < best[0]:
                best = (v, g)
        a = best[1].fuse(cyl(B_JRN[0] / 2, Z["plate_r"] + PLATE_T - Z["gears"] - P.gear_w, Z["gears"] + P.gear_w, cx, cy))
        a = a.fuse(Pos(cx, cy, Z["plate_r"] + PLATE_T) * extrude(RegularPolygon(hexr, 6), amount=GAP_F + 8.0 - 0.3))
        a = a - cyl(R_CORE, 60.0, Z["gears"] - 1.0, cx, cy)
        z1 = Z["plate_r"] + PLATE_T
        b = cyl(B_ECC[0] / 2, Z["disc2"] - z1, z1, cx + E, cy)                              # eccentric 1 (incl. the gap)
        b = b.fuse(cyl(B_ECC[0] / 2, Z["plate_f"] - Z["disc2"], Z["disc2"], cx - E, cy))   # eccentric 2
        b = b.fuse(cyl(B_JRN[0] / 2, PLATE_T, Z["plate_f"], cx, cy))                       # front journal
        b = b - Pos(cx, cy, z1 - 0.1) * extrude(RegularPolygon(hexr + 0.15, 6), amount=GAP_F + 8.0 + 0.1)
        b = b - cyl(R_CORE, 60.0, Z["gears"] - 1.0, cx, cy)
        pieces_a.append(a)
        pieces_b.append(b)
        cores.append(cyl(R_CORE - 0.01, 45.0, Z["gears"] + 2.0, cx, cy))
    p["crank_a"], p["crank_b"], p["crank_core"] = pieces_a, pieces_b, cores

    # input cover (printed): end plate with the NEMA 17 boss and bolt holes, wall round the planets
    r_in = P.r_crank + P.module * P.z_planet / 2 + P.module + 1.5
    cov = ring(r_in, r_out, Z["cover"], Z["plate_r"] - Z["cover"])          # its end face retains the rear main bearing
    end = ring(NEMA["boss"] / 2 + 0.3, r_out, Z["cover"], 4.0)
    s2 = NEMA["bolt_spacing"] / 2
    end = end.cut(*[cyl(1.7, 4.4, Z["cover"] - 0.2, sx * s2, sy * s2) for sx in (-1, 1) for sy in (-1, 1)])
    p["input_cover"] = [cov.fuse(end)]
    motor = Pos(0, 0, -NEMA["body"]) * Box(NEMA["square"], NEMA["square"], NEMA["body"], align=(Align.CENTER, Align.CENTER, Align.MIN))
    motor = motor.fuse(cyl(NEMA["boss"] / 2, 2.0, 0.0), cyl(shaft_d / 2, Z["gears"] + P.gear_w + 1.0, 0.0))
    shaft_flat = Pos(0, shaft_d / 2 - NEMA["flat"] + 1.5, Z["gears"] - 2.0) * Box(shaft_d + 1, 3.0, P.gear_w + 6.0, align=(Align.CENTER, Align.CENTER, Align.MIN))
    p["motor"] = [motor - shaft_flat]
    return p, Z


CHECKS = [("disc", "pin"), ("disc", "ecc_bearing"), ("crank_b", "ecc_bearing"), ("disc", "disc"), ("disc", "pillar"),
          ("disc", "carrier_rear"), ("disc", "carrier_front"), ("disc", "crank_b"), ("pillar", "carrier_rear"),
          ("journal_bearing", "carrier_rear"), ("journal_bearing", "carrier_front"), ("crank_a", "journal_bearing"),
          ("crank_b", "journal_bearing"), ("crank_a", "crank_b"), ("crank_a", "crank_core"), ("crank_b", "crank_core"),
          ("main_bearing", "housing"), ("main_bearing", "carrier_rear"), ("main_bearing", "carrier_front"),
          ("housing", "pin"), ("housing", "disc"), ("sun", "crank_a"), ("crank_a", "crank_a"), ("crank_a", "input_cover"),
          ("crank_a", "carrier_rear"), ("crank_b", "carrier_front"), ("crank_b", "carrier_rear"), ("sun", "motor"),
          ("input_cover", "motor"), ("input_cover", "housing"), ("pin", "carrier_rear"), ("pin", "carrier_front"),
          ("pin", "main_bearing"), ("ecc_bearing", "ecc_bearing")]

PRINTED = ["housing", "disc", "carrier_rear", "carrier_front", "crank_a", "crank_b", "sun", "input_cover"]


def main():
    import cad as steel_cad
    from build123d import Compound, export_step, export_stl
    OUT.mkdir(parents=True, exist_ok=True)
    a = analyse()
    print(f"printable RV: {P.pins} pins Ø{2 * P.pin_r:.0f} on R{P.r_pins:.0f}, {P.lobes} lobes, e {P.e}, K1 {a['k1']:.3f}; "
          f"planetary {P.z_sun}/{P.z_planet} m{P.module} on r{P.r_crank}; ratio {a['ratio']:.1f}:1; "
          f"profile {'no loops' if a['profile_ok'] else 'LOOPS'} (tightest radius {a['tightest_radius']:.2f} mm)")
    print(f"NEMA 17 at {T_MOTOR} N m (assumed) x {a['ratio']:.0f} x eta {ETA} = {a['t_nominal']:.1f} N m nominal output")
    print(f"torque at which each printed feature reaches its allowable ({MAT['name']}, assumed allowables):")
    for k, v in sorted(a["limits"].items(), key=lambda kv: kv[1]):
        val = f"{v:6.1f} N m" if v < 499 else "  > 500 N m (not limiting)"
        print(f"  {k:32} {val}{'   <- governs' if k == a['governing'] else ''}")
    print(f"  lobe contact at the nominal torque: {a['lobe_p0_nominal']:.0f} MPa")
    rows, cp = ring_options()
    print("\nring pins, efficiency (cycloidal/efficiency.py's contact model; ranking robust, absolute values uncalibrated):")
    for r in rows:
        print(f"  {r['name']:42} mu {r['mu']:5.3f}  cycloid stage {r['eta_cyclo']:.2f}  whole reducer {r['eta']:.2f}  "
              f"back-driving {r['eta_back']:.2f}  output from 0.4 N m in: {r['t_out']:4.1f} N m")
    print(f"  a printed tooth is as soft as the lobe: E* {cp['e_star_printed_tooth']:.0f} vs {cp['e_star_steel_pin']:.0f} MPa, "
          f"so contact pressure x{cp['pressure_ratio']:.2f} at the same load")
    print(f"bearings at the nominal torque (equivalent crank-bearing load {a['bearing_force_eq']:.0f} N, crank {a['crank_rpm']:.0f} rpm): "
          + ", ".join(f"{k} L10 {v:.3g} h" for k, v in a["lives_h"].items()))

    p, Z = build()
    print(f"\nassembly: Ø{2 * (B_MAIN[1] / 2 + 5):.0f} mm, {Z['end'] + 2 + NEMA['body']:.0f} mm long with the motor")
    bad = 0
    print("interference, mm3:")
    for x, y in CHECKS:
        va, vb = p[x], p[y]
        pairs = ([(va[i], va[j]) for i in range(len(va)) for j in range(i + 1, len(va))] if x == y
                 else [(u, w) for u in va for w in vb])
        v = sum(steel_cad.shared(u, w) for u, w in pairs)
        bad += v >= 1e-3
        if v >= 1e-3:
            print(f"  {x:>16} / {y:<16} {v:10.4f}  <--")
    print(f"  {len(CHECKS)} pairs checked, {bad} interfering")
    failed = steel_cad.check_parts(p)
    for name, solids in p.items():
        export_stl(Compound(children=solids), str(OUT / f"{name}.stl"), tolerance=0.02, angular_tolerance=0.2)
    # print files: one of each printed part, as placed (flat faces down in the slicer)
    (OUT / "print").mkdir(exist_ok=True)
    for name in PRINTED:
        export_stl(p[name][0], str(OUT / "print" / f"{name}.stl"), tolerance=0.01, angular_tolerance=0.1)
    export_step(steel_cad.named_assembly(p), str(OUT / "assembly.step"))
    bom = dict(
        bought={"6802-2RS ball bearing 15x24x5": 12, "6800-2RS ball bearing 10x19x5": 6, "6814-2RS ball bearing 70x90x10": 2,
                "DIN 6325 hardened dowel pin 5x24 (pin ring)": P.pins, "DIN 6325 hardened dowel pin 5x45 (pillars)": 3,
                "DIN 6325 hardened dowel pin 5x45 (crankshaft cores)": 3, "NEMA 17 stepper, 5 mm D shaft": 1,
                "M3 x 8 screws (motor)": 4},
        printed={"housing": 1, "disc (print 2; flip one over at assembly)": 2, "carrier_rear": 1, "carrier_front": 1,
                 "crank_a (planet + rear journal)": 3, "crank_b (eccentrics + front journal)": 3, "sun": 1, "input_cover": 1})
    (OUT / "bom.json").write_text(json.dumps(dict(analysis=a, bom=bom, ring_options=rows, contact=cp), indent=1, default=float))
    print(f"\n{len(failed)} failed solids. wrote {OUT}/assembly.step, *.stl, print/*.stl, bom.json")
    if bad or failed:
        for f in failed:
            print("  ", f)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
