"""
CAD of the RV actuator (design.py's numbers, out/design.json), as an assembly with
interference and per-solid checks.

Axis z, input (motor) end at z = 0, output flange at the far end. Coaxial, so parts are
built in place from the stations below.

  motor           stator ring with teeth, rotor (yoke + magnets) on the input shaft
  input_shaft     motor rotor to the sun gear, on two bearings
  sun, planet x3  first stage; each planet sits on a crankshaft end
  crankshaft x3   planet seat, rear journal, two eccentrics 180 deg apart, front journal
  roller x84      needle rollers, 14 per eccentric, running on the eccentric and in the disc bore
  disc x2         cycloid discs, offset +e and -e
  pin x40         in half-round grooves in the pin case
  carrier_front   output flange with three pillars through the discs
  carrier_rear    bolts to the pillars; both flanges carry the crank journal bearings
  housing         pin case + main bearing seats; gear_cover and motor can at the input end
  bearings        main (x2), crank journal (x6), input shaft (x2): plain rings, stand-ins

Stand-ins, stated: bearings are solid rings; gear teeth are build123d involutes, phased by
trial against the interference check; no fasteners, seals, oil or encoder detail.

    python rv/cad.py          (or: make rv-cad)
"""

import json
import sys
from math import cos, radians, sin, tau
from pathlib import Path

from build123d import Align, Compound, Cylinder, Pos, Rot, export_step, export_stl
from build123d import BuildLine, BuildPart, BuildSketch, Solid, Spline, extrude, make_face

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402

OUT = HERE / "out" / "cad"
J = json.loads((HERE / "out" / "design.json").read_text())
RV = dz.RV(**J["rv"])
MOT = J["motor"]

E = RV.e
R_ECC, R_ROLL, R_CRANK = RV.ecc_r, RV.roller_r, RV.r_crank
BRG = J["bearings"]
ECC_B = BRG["eccentric"]                           # RNA4906: outer ring pressed into the disc, needles on the eccentric
R_BORE = ECC_B["D"] / 2                            # disc bore = the bearing's OD
R_RING_IN = R_ECC + 2 * R_ROLL                     # outer ring raceway (needles between it and the eccentric)
L_NEEDLE = RV.roller_len or RV.disc_t
B_JOURNAL = (BRG["journal"]["d"] / 2, BRG["journal"]["D"] / 2, BRG["journal"]["B"])     # 32005 X
R_JOURNAL = B_JOURNAL[0]
R_PILLAR = 9.0
R_PILLAR_HOLE = R_PILLAR + E + 0.5
B_MAIN = (BRG["main"]["d"] / 2, BRG["main"]["D"] / 2, BRG["main"]["B"])               # 61830, one per flange
R_FLANGE = B_MAIN[0]                                 # carrier flange radius = main bearing bore
R_PIN_SEAT = RV.r_pins
R_HOUSING_IN = RV.r_pins                            # pins half-sunk in the pin case
R_HOUSING = B_MAIN[1] + 9.0
DISC_T, GAP = RV.disc_t, 1.0
PIN_L = J["pin"]["l"]                               # DIN 6325 dowel pins, 8 x 40
GAP_F = (PIN_L - 2 * DISC_T - GAP) / 2              # flange to disc: the pins exactly span both discs and gaps
MOD = 0.03                                          # mm, equidistant modification of the disc profile (30 um)
FLANGE_T = max(B_JOURNAL[2] + 3.0, B_MAIN[2])

# planetary
M_GEAR, Z_SUN, Z_PL, GEAR_W = RV.module, RV.z_sun, RV.z_planet, RV.gear_w
R_PL_TIP = M_GEAR * Z_PL / 2 + M_GEAR
R_SHAFT = 10.0
B_INPUT = (BRG["input"]["d"] / 2, BRG["input"]["D"] / 2, BRG["input"]["B"])         # 6004

# motor (from design.json)
M_OD, M_STACK, M_BORE, M_SLOT = MOT["od"], MOT["stack"], MOT["bore"], MOT["slot_depth"]
M_OVER = M_SLOT / 2 + 2

# ---- axial stations, input end first ---------------------------------------------------------
Z = {}
Z["motor_cap"] = 0.0
Z["bearing_in_rear"] = 5.0
Z["stack"] = Z["bearing_in_rear"] + B_INPUT[2] + 3 + M_OVER
Z["bearing_in_front"] = Z["stack"] + M_STACK + M_OVER + 3
Z["gears"] = Z["bearing_in_front"] + B_INPUT[2] + 4
Z["flange_rear"] = Z["gears"] + GEAR_W + 4
Z["disc1"] = Z["flange_rear"] + FLANGE_T + GAP_F
Z["disc2"] = Z["disc1"] + DISC_T + GAP
Z["flange_front"] = Z["disc2"] + DISC_T + GAP_F
Z["end"] = Z["flange_front"] + FLANGE_T


def ring(r_in, r_out, z0, h, x=0.0, y=0.0):
    c = Cylinder(r_out, h, align=(Align.CENTER, Align.CENTER, Align.MIN))
    if r_in > 0:
        c -= Cylinder(r_in, h + 0.2, align=(Align.CENTER, Align.CENTER, Align.MIN)).moved(Pos(0, 0, -0.1))
    return Pos(x, y, z0) * c


# Cranks at 90, 210, 330 deg: one crank on the axis at right angles to the profile's symmetry
# line (a pin valley sits on +x). Then the hole pattern is symmetric about the y axis while
# the 39-lobe profile, mirrored about y, shifts half a lobe, which is exactly what the
# 180-deg eccentric needs. So disc 2 is disc 1 flipped over: one part number. (With the
# cranks at 0/120/240 no rotation does both: 120 deg is never a whole number of lobes plus a half.)
CRANK0 = 90.0


def crank_xy(k):
    a = radians(CRANK0) + tau * k / RV.n_crank
    return R_CRANK * cos(a), R_CRANK * sin(a)


def pillar_xy(k):
    a = radians(CRANK0) + tau * (k + 0.5) / RV.n_crank
    return R_CRANK * cos(a), R_CRANK * sin(a)


# ---- cycloid stage ---------------------------------------------------------------------------------
def disc_local():
    """disc in its own frame (centre at the origin)"""
    # One periodic B-spline through 64 points a lobe: within 0.8 um of the true profile, and
    # a STEP outline of one curve rather than thousands of line segments (which made the file
    # 37 MB and too slow for viewers such as CAD Assistant).
    x, y = dz.disc_profile(RV, steps=RV.lobes * 64, mod=MOD)   # equidistant profile modification: running clearance
    pts = list(zip(x.tolist(), y.tolist()))
    with BuildPart() as p:
        with BuildSketch():
            with BuildLine():
                Spline(*pts, periodic=True)
            make_face()
        extrude(amount=DISC_T)
    d = p.part
    holes = [Pos(*crank_xy(k), -0.1) * Cylinder(R_BORE, DISC_T + 0.2, align=(Align.CENTER, Align.CENTER, Align.MIN))
             for k in range(RV.n_crank)]
    holes += [Pos(*pillar_xy(k), -0.1) * Cylinder(R_PILLAR_HOLE, DISC_T + 0.2, align=(Align.CENTER, Align.CENTER, Align.MIN))
              for k in range(RV.n_crank)]
    holes.append(Pos(0, 0, -0.1) * Cylinder(14.0, DISC_T + 0.2, align=(Align.CENTER, Align.CENTER, Align.MIN)))
    return d.cut(*holes)


def discs():
    from build123d import Plane, mirror
    d = disc_local()
    d1 = Pos(E, 0, Z["disc1"]) * d
    d2 = Pos(-E, 0, Z["disc2"]) * mirror(d, about=Plane.YZ)       # the same disc, flipped over
    return [d1, d2]


def pins():
    return [Pos(R_PIN_SEAT * cos(tau * k / RV.pins), R_PIN_SEAT * sin(tau * k / RV.pins), Z["disc1"] - GAP_F)
            * Cylinder(RV.pin_r, PIN_L, align=(Align.CENTER, Align.CENTER, Align.MIN))
            for k in range(RV.pins)]


def crankshaft(k):
    cx, cy = crank_xy(k)
    parts = [ring(0, 10.0, Z["gears"] - 1, GEAR_W + 2, cx, cy),                       # planet seat (spline stand-in)
             ring(0, R_JOURNAL, Z["gears"] + GEAR_W, Z["disc1"] - Z["gears"] - GEAR_W + 0.5, cx, cy),
             ring(0, R_ECC, Z["disc1"] - GAP / 2, DISC_T + GAP / 2 + 0.5, cx + E, cy),        # eccentric 1 (+e)
             ring(0, R_ECC, Z["disc2"] - 0.5, DISC_T + GAP / 2 + 0.5, cx - E, cy),            # eccentric 2 (-e)
             ring(0, R_JOURNAL, Z["disc2"] + DISC_T - 0.5, Z["end"] - Z["disc2"] - DISC_T - 1.5, cx, cy)]
    return parts[0].fuse(*parts[1:])


def rollers():
    """the RNA4906 needles, running directly on the hardened eccentric (no inner ring)"""
    out = []
    rp = R_ECC + R_ROLL
    for k in range(RV.n_crank):
        cx, cy = crank_xy(k)
        for side, z0 in ((1, Z["disc1"]), (-1, Z["disc2"])):
            for i in range(RV.n_rollers):
                a = tau * i / RV.n_rollers
                out.append(Pos(cx + side * E + rp * cos(a), cy + rp * sin(a), z0 + (DISC_T - L_NEEDLE) / 2)
                           * Cylinder(R_ROLL - 0.005, L_NEEDLE, align=(Align.CENTER, Align.CENTER, Align.MIN)))
    return out


def ecc_rings():
    """the RNA4906 outer rings, pressed into the disc bores"""
    out = []
    for k in range(RV.n_crank):
        cx, cy = crank_xy(k)
        for side, z0 in ((1, Z["disc1"]), (-1, Z["disc2"])):
            out.append(ring(R_RING_IN, R_BORE, z0, DISC_T, cx + side * E, cy))
    return out


def carrier_front():
    f = ring(15.0, R_FLANGE, Z["flange_front"], FLANGE_T)
    for k in range(RV.n_crank):
        f -= Pos(*crank_xy(k), Z["flange_front"] - 0.1) * Cylinder(B_JOURNAL[1], FLANGE_T + 0.2,
                                                                       align=(Align.CENTER, Align.CENTER, Align.MIN))
    pillars = [Pos(*pillar_xy(k), Z["flange_rear"] + FLANGE_T) * Cylinder(R_PILLAR, Z["flange_front"] - Z["flange_rear"] - FLANGE_T + 1.0,
                                                                         align=(Align.CENTER, Align.CENTER, Align.MIN))
               for k in range(RV.n_crank)]
    return f.fuse(*pillars)


def carrier_rear():
    f = ring(15.0, R_FLANGE, Z["flange_rear"], FLANGE_T)
    for k in range(RV.n_crank):
        f -= Pos(*crank_xy(k), Z["flange_rear"] - 0.1) * Cylinder(B_JOURNAL[1], FLANGE_T + 0.2,
                                                                      align=(Align.CENTER, Align.CENTER, Align.MIN))
    return f


def journal_bearings():
    out = []
    for k in range(RV.n_crank):
        cx, cy = crank_xy(k)
        for zf in (Z["flange_rear"] + (FLANGE_T - B_JOURNAL[2]) / 2, Z["flange_front"] + (FLANGE_T - B_JOURNAL[2]) / 2):
            out.append(ring(B_JOURNAL[0], B_JOURNAL[1], zf, B_JOURNAL[2], cx, cy))
    return out


def main_bearings():
    return [ring(B_MAIN[0], B_MAIN[1], Z["flange_rear"] + FLANGE_T - B_MAIN[2], B_MAIN[2]),
            ring(B_MAIN[0], B_MAIN[1], Z["flange_front"], B_MAIN[2])]


def housing():
    z0, z1 = Z["flange_rear"] + FLANGE_T - B_MAIN[2], Z["flange_front"] + B_MAIN[2]
    h = ring(B_MAIN[1], R_HOUSING, z0, z1 - z0)
    h = h.fuse(ring(R_HOUSING_IN, R_HOUSING - 1, Z["disc1"] - GAP_F, PIN_L))   # pin case
    grooves = [Pos(R_PIN_SEAT * cos(tau * k / RV.pins), R_PIN_SEAT * sin(tau * k / RV.pins), Z["disc1"] - GAP_F - 0.1)
               * Cylinder(RV.pin_r, PIN_L + 0.2, align=(Align.CENTER, Align.CENTER, Align.MIN))
               for k in range(RV.pins)]
    return h.cut(*grooves)


def gear_cover():
    z0 = Z["bearing_in_front"]
    r_in = R_PL_TIP + R_CRANK + 2.0
    c = ring(r_in, r_in + 4.0, z0, Z["flange_rear"] + FLANGE_T - B_MAIN[2] - z0)
    wall = ring(B_INPUT[1], r_in + 4.0, z0, 4.0)
    return c.fuse(wall)


def motor_parts():
    r_so = M_OD / 2 - 2.0
    can = ring(r_so, M_OD / 2, Z["motor_cap"], Z["bearing_in_front"] - Z["motor_cap"])
    cap = ring(B_INPUT[1], M_OD / 2, Z["motor_cap"], 5.0)
    r_si = M_BORE / 2 + 2.5 + 2.5 + 0.5           # bore + rotor yoke + magnet + gap (rollerscrew Motor defaults)
    stator = ring(r_si, r_so, Z["stack"], M_STACK)
    coils = [ring(r_si + 1.5, r_so - 3.0, Z["stack"] - M_OVER + 0.5, M_OVER - 0.5),           # end-winding envelopes
             ring(r_si + 1.5, r_so - 3.0, Z["stack"] + M_STACK, M_OVER - 0.5)]
    rotor = ring(R_SHAFT, r_si - 0.5, Z["stack"], M_STACK)
    return dict(motor_can=[can.fuse(cap)], stator=[stator], end_windings=coils, rotor=[rotor])


def input_parts():
    shaft = ring(0, R_SHAFT, Z["bearing_in_rear"], Z["gears"] - Z["bearing_in_rear"])
    brgs = [ring(B_INPUT[0], B_INPUT[1], Z["bearing_in_rear"], B_INPUT[2]),
            ring(B_INPUT[0], B_INPUT[1], Z["bearing_in_front"], B_INPUT[2])]
    return shaft, brgs


SUN_TURN, PLANET_TURN = 0.0, [0.0, 0.0, 0.0]
BACKLASH_FLANK = 0.04                               # mm per flank, taken off the sun


def spur_gear(module, teeth, thickness, bore, pressure_angle=20.0):
    """involute spur gear, base at z = 0 (bd_warehouse). Above ~41 teeth at 20 deg the standard
    1.25 m dedendum puts the root outside the base circle, which the library refuses; the
    dedendum is deepened just enough to bring the root inside (more root clearance, same flank)."""
    from math import cos as _cos
    from bd_warehouse.gear import SpurGear
    r = module * teeth / 2
    ded = max(1.25 * module, r * (1 - _cos(radians(pressure_angle))) + 0.3)     # +0.3 leaves room for the root fillet
    g = SpurGear(module=module, tooth_count=teeth, thickness=thickness, pressure_angle=pressure_angle,
                 root_fillet=0.2 * module, dedendum=ded)
    if bore > 0:
        g = g - Cylinder(bore / 2, 2 * thickness)
    return Pos(0, 0, thickness / 2) * g


def sun():
    """pinion cut on the end of the input shaft (a separate solid here, butting the shaft).
    Every face offset 0.04 mm inward: backlash, since standard gears at the standard centre
    distance have none (and the flanks would touch)."""
    from build123d import Kind, offset
    from OCP.BRepBuilderAPI import BRepBuilderAPI_NurbsConvert
    g = offset(spur_gear(M_GEAR, Z_SUN, GEAR_W, bore=0), amount=-BACKLASH_FLANK, kind=Kind.INTERSECTION)
    # the offset leaves OFFSET_SURFACE entities, which many STEP readers handle badly: convert to NURBS
    g = Solid(BRepBuilderAPI_NurbsConvert(g.solids()[0].wrapped, True).Shape())
    return Pos(0, 0, Z["gears"]) * Rot(0, 0, SUN_TURN) * g


def planet(k, turn=None):
    t = PLANET_TURN[k] if turn is None else turn
    return Pos(*crank_xy(k), Z["gears"]) * Rot(0, 0, t) * spur_gear(M_GEAR, Z_PL, GEAR_W, bore=20.0)


def assembly():
    shaft, brg_in = input_parts()
    p = dict(housing=[housing()], pin=pins(), disc=discs(), crankshaft=[crankshaft(k) for k in range(RV.n_crank)],
             roller=rollers(), ecc_ring=ecc_rings(), carrier_front=[carrier_front()], carrier_rear=[carrier_rear()],
             journal_bearing=journal_bearings(), main_bearing=main_bearings(), sun=[sun()],
             planet=[planet(k) for k in range(RV.n_crank)], input_shaft=[shaft], input_bearing=brg_in,
             gear_cover=[gear_cover()])
    p.update(motor_parts())
    return p


def shared(a, b):
    ba, bb = a.bounding_box(), b.bounding_box()
    if any(tuple(ba.min)[i] > tuple(bb.max)[i] or tuple(bb.min)[i] > tuple(ba.max)[i] for i in range(3)):
        return 0.0
    r = a.intersect(b)
    if r is None:
        return 0.0
    return sum(x.volume for x in r) if isinstance(r, (list, tuple)) or hasattr(r, "__iter__") and not hasattr(r, "volume") else r.volume


def phase_by_trial():
    """each planet's turn, chosen by least interference with the sun (coarse, then fine)"""
    global PLANET_TURN
    s = sun()
    pitch = 360.0 / Z_PL
    turns = []
    for k in range(RV.n_crank):
        v, t = min((shared(s, planet(k, t)), t) for t in [pitch * i / 12 for i in range(12)])
        v, t = min((shared(s, planet(k, tt)), tt) for tt in [t + pitch * j / 96 for j in range(-6, 7)])
        turns.append(t)
        print(f"planet {k + 1}: turned {t:.3f} deg, {v:.4f} mm3 against the sun")
    PLANET_TURN = turns


CHECKS = [("disc", "pin"), ("disc", "roller"), ("crankshaft", "roller"), ("roller", "roller"), ("ecc_ring", "roller"),
          ("ecc_ring", "disc"), ("ecc_ring", "crankshaft"), ("pin", "main_bearing"), ("pin", "carrier_front"),
          ("pin", "carrier_rear"), ("disc", "carrier_front"),
          ("disc", "carrier_rear"), ("disc", "disc"), ("disc", "crankshaft"), ("crankshaft", "journal_bearing"),
          ("carrier_front", "journal_bearing"), ("carrier_rear", "journal_bearing"), ("main_bearing", "housing"),
          ("main_bearing", "carrier_front"), ("main_bearing", "carrier_rear"), ("housing", "pin"), ("housing", "disc"),
          ("sun", "planet"), ("planet", "planet"), ("planet", "gear_cover"), ("planet", "carrier_rear"),
          ("crankshaft", "carrier_front"), ("crankshaft", "carrier_rear"), ("carrier_front", "carrier_rear"),
          ("input_shaft", "sun"), ("rotor", "stator"), ("stator", "motor_can"), ("input_bearing", "input_shaft"),
          ("gear_cover", "housing"), ("crankshaft", "planet")]


def interference(p):
    print("\ninterference, mm3 shared (press fits and seats read as face contact, 0):")
    bad = 0
    for a, b in CHECKS:
        va, vb = p[a], p[b]
        if a == b:
            pairs = [(va[i], va[j]) for i in range(len(va)) for j in range(i + 1, len(va))]
            if a == "roller":
                pairs = [(va[i], va[i + 1]) for i in range(RV.n_rollers - 1)]
        else:
            pairs = [(x, y) for x in va for y in vb]
        v = sum(shared(x, y) for x, y in pairs)
        flag = "" if v < 1e-3 else "  <--"
        bad += v >= 1e-3
        print(f"  {a:>16} / {b:<16} {v:10.4f}{flag}")
    return bad


def check_parts(p):
    import tempfile

    import trimesh
    bad = []
    print("\npart checks (each solid: watertight, one body):")
    with tempfile.TemporaryDirectory() as tmp:
        for name, solids in p.items():
            ok = 0
            for k, s in enumerate(solids):
                f = f"{tmp}/{name}{k}.stl"
                export_stl(s, f, tolerance=0.02, angular_tolerance=0.2)
                m = trimesh.load(f)
                nb = len(m.split(only_watertight=False))
                if m.is_watertight and nb == 1:
                    ok += 1
                else:
                    bad.append(f"{name}[{k}]: watertight={m.is_watertight} bodies={nb}")
            bb = Compound(children=solids).bounding_box()
            print(f"  {name:16} {ok:3d}/{len(solids):<3} ok  bbox {bb.size.X:6.1f} x {bb.size.Y:6.1f} x {bb.size.Z:6.1f}")
    return bad


def named_assembly(p):
    """a STEP assembly: rv_actuator / <part type> / <part type>_<n>, so viewers show a tree"""
    groups = []
    for name, solids in p.items():
        kids = []
        for k, s in enumerate(solids):
            s = s.solids()[0] if hasattr(s, "solids") and len(s.solids()) == 1 else s
            s.label = f"{name}_{k + 1}" if len(solids) > 1 else name
            kids.append(s)
        groups.append(Compound(children=kids, label=name))
    return Compound(children=groups, label="rv_actuator")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"RV {RV.ratio:.2f}:1, Ø{2 * R_HOUSING:.0f} x {Z['end']:.0f} mm overall; stations: "
          + ", ".join(f"{k} {v:.1f}" for k, v in Z.items()))
    phase_by_trial()
    p = assembly()
    for name, solids in p.items():
        export_stl(Compound(children=solids), str(OUT / f"{name}.stl"), tolerance=0.02, angular_tolerance=0.2)
    export_step(named_assembly(p), str(OUT / "assembly.step"))
    n_bad = interference(p)
    bad = check_parts(p)
    if bad:
        print("FAILED part checks:\n  " + "\n  ".join(bad))
    print(f"\n{n_bad} interfering pair types; {len(bad)} failed solids. wrote {OUT}/*.stl, assembly.step")
    if n_bad or bad:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
