"""
CAD of the shallow-slot liquid-cooled motor (design.SHALLOW), as an assembly.

Every dimension comes from design.SHALLOW or from the constants below; nothing is
re-typed in a part. The motor is coaxial, so the mates are axial stations: each part
is built with its seat at z = 0 and placed at a station from STATIONS.

Parts
  stator_core      24 teeth, parallel-sided slots, tooth tips
  coil             one tooth's four turns of hollow rectangular conductor (x24)
  jacket_inner     sleeve the stator presses into; lands at each end close the channel
  jacket_outer     outer skin; with jacket_inner it forms the annular coolant channel
  rotor_yoke, magnet (x22), rotor_carrier
  manifold         insulating ring with supply and return plenums, one port pair per coil
  drive_plate      ring cold plate carrying the drive, first on the coolant loop
  cap_rear, cap_front, shaft
  bearing (x2)     stock 6903, stand-in

Stand-ins, stated: each coil turn is a closed racetrack (no turn-to-turn crossover or
leads); the bearing and shaft arrangement is a placeholder with no gear stage; no
fasteners, seals or bus bars.

    python qdd-liquid/cad.py        (or: make qdd-liquid-cad)
"""

import sys
from math import asin, atan2, cos, degrees, hypot, pi, radians, sin, tau
from pathlib import Path

from build123d import (Align, Axis, Box, BuildLine, BuildPart, BuildSketch, Circle, Compound,
                       Cylinder, Line, Location, Locations, Mode, Plane, PolarLocations, Polygon,
                       Pos, RadiusArc, Rectangle, Rot, ThreePointArc, add, export_step, export_stl,
                       extrude, make_face, sweep)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from design import SHALLOW as D  # noqa: E402

OUT = HERE / "out" / "cad"

# ---- radii, mm: all from the design ---------------------------------------------
R_GAP = D.r_gap * 1e3                       # magnet outer radius
R_BORE = R_GAP + D.air_gap                  # stator bore
R_SLOT_A, R_SLOT_B = (x * 1e3 for x in D.slot_radii)
R_STATOR = R_SLOT_B + D.stator_yoke
R_MAG_A = R_GAP - D.magnet
R_ROTOR_IN = R_MAG_A - D.rotor_yoke
R_OUT = D.od / 2
STACK = D.stack
SLOT_PITCH, POLE_PITCH = tau / D.slots, tau / D.poles
OPEN_FRAC, MAGNET_ARC = 0.2, 0.85           # as solved in fea.Motor

# ---- jacket: inner wall, channel, outer wall inside D.jacket ----------------------
J_WALL = 1.0
J_CHANNEL = D.jacket - 2 * J_WALL
LAND = 5.0                                  # sealing land at each end of the channel

# ---- coil -----------------------------------------------------------------------
CU_W, CU_H = D.cond_w - 2 * D.insulation, D.cond_h - 2 * D.insulation
LEG_V = D.slot_w / 2 - D.liner - D.cond_w / 2           # leg centre, off the slot axis
COIL_EXT = 1.0                              # straight run past the stack before the bend


def turn_radius(i):
    """Radius of turn i's leg centres (0 = nearest the bore)."""
    return R_SLOT_A + D.liner + D.cond_h * (i + 0.5)


def leg_xy(i, side):
    """Leg centre of turn i for the coil on the tooth at angle 0. side -1 is the leg in
    the slot below the tooth (its upper half), +1 the leg in the slot above."""
    a = side * SLOT_PITCH / 2               # that slot's axis
    u, v = turn_radius(i), -side * LEG_V
    return u * cos(a) - v * sin(a), u * sin(a) + v * cos(a)


def bend_radius(i):
    (x0, y0), (x1, y1) = leg_xy(i, -1), leg_xy(i, 1)
    return hypot(x1 - x0, y1 - y0) / 2


OVERHANG = COIL_EXT + bend_radius(D.turns - 1) + CU_W / 2     # end winding past the stack

# ---- axial stations (z of each part's seat), rear to front -----------------------
CAP = 4.0
DRIVE_ZONE = 8.0                            # drive plate and board
MANIFOLD_H = 8.0
END_CLEAR = 1.5
SHAFT_R, BRG_OD, BRG_W = 8.5, 30.0, 7.0     # 6903: 17 x 30 x 7
WEB = 4.0

STATIONS = {}
STATIONS["cap_rear"] = 0.0
STATIONS["drive_plate"] = CAP
STATIONS["manifold"] = CAP + DRIVE_ZONE
STATIONS["stack"] = STATIONS["manifold"] + MANIFOLD_H + END_CLEAR + OVERHANG
STATIONS["cap_front"] = STATIONS["stack"] + STACK + OVERHANG + END_CLEAR
LENGTH = STATIONS["cap_front"] + CAP + 1.0  # front cap is CAP thick plus a bearing boss
STATIONS["jacket"] = CAP


def ring(r_in, r_out, h):
    return Cylinder(r_out, h, align=(Align.CENTER, Align.CENTER, Align.MIN)) \
        - Cylinder(r_in, h, align=(Align.CENTER, Align.CENTER, Align.MIN))


# ---- parts (seat at z = 0) -------------------------------------------------------
def stator_core():
    with BuildPart() as p:
        with BuildSketch():
            Circle(R_STATOR)
            Circle(R_BORE, mode=Mode.SUBTRACT)
            with PolarLocations(0, D.slots):
                with Locations(Pos((R_SLOT_A + R_SLOT_B) / 2, 0)):
                    Rectangle(D.slot_depth, D.slot_w, mode=Mode.SUBTRACT)
                opening = 2 * R_BORE * sin(OPEN_FRAC * SLOT_PITCH / 2)
                with Locations(Pos((R_BORE + R_SLOT_A) / 2, 0)):
                    Rectangle(D.tooth_tip + 1.0, opening, mode=Mode.SUBTRACT)
        extrude(amount=STACK)
    return p.part


def coil_turn(i):
    """One closed racetrack turn around the tooth at angle 0, z = 0 at the stack's rear
    face. The two legs lie square in their own slots; the end bends are semicircles in
    the plane through both legs, so the conductor twists half a slot pitch in each bend."""
    (x0, y0), (x1, y1) = leg_xy(i, -1), leg_xy(i, 1)
    r = bend_radius(i)
    cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
    ang = atan2(y1 - y0, x1 - x0)
    length = STACK + 2 * COIL_EXT
    legs = []
    for side, (x, y) in ((-1, (x0, y0)), (1, (x1, y1))):
        bar = Box(CU_H, CU_W, length) - Cylinder(D.cond_bore / 2, length)     # radial x across x axial
        legs.append(Pos(x, y, STACK / 2) * Rot(0, 0, degrees(side * SLOT_PITCH / 2)) * bar)
    # bends: local X along the chord between the legs, local Y along the motor axis
    plane = Plane(origin=(cx, cy, STACK / 2), x_dir=(cos(ang), sin(ang), 0), z_dir=(sin(ang), -cos(ang), 0))
    half = length / 2
    bends = []
    for sgn in (1, -1):
        with BuildPart() as p:
            with BuildLine(plane) as path:
                ThreePointArc((-r, sgn * half), (0, sgn * (half + r)), (r, sgn * half))
            with BuildSketch(Plane(origin=path.wire() @ 0, z_dir=path.wire() % 0, x_dir=plane.x_dir)):
                Rectangle(CU_W, CU_H)
                Circle(D.cond_bore / 2, mode=Mode.SUBTRACT)
            sweep(path=path.wire())
        bends.append(p.part)
    return legs[0] + legs[1] + bends[0] + bends[1]


def coil():
    return Compound(children=[coil_turn(i) for i in range(D.turns)])


def rotor_yoke():
    return ring(R_ROTOR_IN, R_MAG_A, STACK)


def magnet():
    """One pole, centred on angle 0."""
    half = MAGNET_ARC * POLE_PITCH / 2
    with BuildPart() as p:
        with BuildSketch():
            with BuildLine():
                RadiusArc((R_GAP * cos(half), -R_GAP * sin(half)), (R_GAP * cos(half), R_GAP * sin(half)), -R_GAP)
                Line((R_GAP * cos(half), R_GAP * sin(half)), (R_MAG_A * cos(half), R_MAG_A * sin(half)))
                RadiusArc((R_MAG_A * cos(half), R_MAG_A * sin(half)), (R_MAG_A * cos(half), -R_MAG_A * sin(half)), R_MAG_A)
                Line((R_MAG_A * cos(half), -R_MAG_A * sin(half)), (R_GAP * cos(half), -R_GAP * sin(half)))
            make_face()
        extrude(amount=STACK)
    return p.part


def rotor_carrier():
    """Web at the front of the rotor yoke down to a hub on the shaft. Seat: yoke front face."""
    web = ring(SHAFT_R, R_MAG_A, WEB)
    hub = Pos(0, 0, -6) * ring(SHAFT_R, SHAFT_R + 4, 6)
    return web + hub


def jacket_inner():
    """Seat: rear end. The stator presses into the bore; the lands seal to jacket_outer."""
    h = LENGTH - 2 * CAP
    body = ring(R_STATOR, R_STATOR + J_WALL, h)
    lands = ring(R_STATOR, R_OUT - J_WALL, LAND) + Pos(0, 0, h - LAND) * ring(R_STATOR, R_OUT - J_WALL, LAND)
    return body + lands


def jacket_outer():
    h = LENGTH - 2 * CAP
    shell = ring(R_OUT - J_WALL, R_OUT, h)
    ports = []
    for z, a in ((LAND + 4, 0), (h - LAND - 4, 180)):          # inlet at the rear, outlet at the front
        boss = Rot(0, 0, a) * (Pos(R_OUT + 3, 0, z) * Rot(0, 90, 0) * ring(2.0, 4.0, 6))
        ports.append(boss)
    bore = [Rot(0, 0, a) * (Pos(R_OUT, 0, z) * Rot(0, 90, 0) * Cylinder(2.0, 6)) for z, a in
            ((LAND + 4, 0), (h - LAND - 4, 180))]
    return shell + ports[0] + ports[1] - bore[0] - bore[1]


def manifold():
    """Insulating ring behind the rear end windings: an outer supply plenum and an inner
    return plenum, each open to every coil through one port. Seat: rear face."""
    r_in, r_out = R_SLOT_A - 2.0, R_SLOT_B + 2.0
    body = ring(r_in, r_out, MANIFOLD_H)
    mid = (r_in + r_out) / 2
    supply = Pos(0, 0, 1.5) * ring(mid + 1.0, r_out - 1.5, MANIFOLD_H - 4.0)
    ret = Pos(0, 0, 1.5) * ring(r_in + 1.5, mid - 1.0, MANIFOLD_H - 4.0)
    body = body - supply - ret
    with BuildPart() as ports:
        with PolarLocations(0, D.slots, start_angle=degrees(SLOT_PITCH / 2)):
            for r in ((mid + 1.0 + r_out - 1.5) / 2, (r_in + 1.5 + mid - 1.0) / 2):
                with Locations(Pos(r, 0, MANIFOLD_H - 1.5)):
                    Cylinder(D.cond_bore / 2 + 0.4, 6)
    return body - ports.part


def drive_plate():
    """Ring cold plate for the drive, with one annular coolant passage. Seat: rear face."""
    body = ring(SHAFT_R + 10, R_STATOR - 1.0, DRIVE_ZONE - 2.0)
    passage = Pos(0, 0, 1.5) * ring(R_SLOT_A, R_SLOT_A + 6.0, DRIVE_ZONE - 5.0)
    return body - passage


def cap(front):
    """End cap with a bearing seat. Seat: its inner face at z = 0, body toward +z."""
    plate = ring(BRG_OD / 2 - 2.0, R_OUT, CAP)
    boss = ring(BRG_OD / 2, BRG_OD / 2 + 3.0, BRG_W) if front else Pos(0, 0, 0) * ring(BRG_OD / 2, BRG_OD / 2 + 3.0, BRG_W)
    return plate + (Pos(0, 0, CAP - BRG_W) * boss if front else boss)


def bearing():
    return ring(SHAFT_R, BRG_OD / 2, BRG_W)


def shaft():
    return ring(SHAFT_R - 4.0, SHAFT_R, LENGTH + 6.0)


# ---- assembly --------------------------------------------------------------------
COLORS = {"stator_core": "#5b6470", "coil": "#c87a2f", "rotor_yoke": "#34404d", "magnet": "#7aa6c2",
          "rotor_carrier": "#8d99a6", "jacket_inner": "#b9c2cc", "jacket_outer": "#d9dee4",
          "manifold": "#d8c9a3", "drive_plate": "#3f7f5f", "cap_rear": "#aab2bb", "cap_front": "#aab2bb",
          "bearing": "#e6e6e6", "shaft": "#777f88"}


def assembly():
    """{name: [placed solids]}."""
    zs = STATIONS["stack"]
    one_coil = coil()
    parts = {
        "stator_core": [Pos(0, 0, zs) * stator_core()],
        "coil": [Pos(0, 0, zs) * Rot(0, 0, degrees((k + 0.5) * SLOT_PITCH)) * one_coil for k in range(D.slots)],
        "rotor_yoke": [Pos(0, 0, zs) * rotor_yoke()],
        "magnet": [Pos(0, 0, zs) * Rot(0, 0, degrees(k * POLE_PITCH)) * magnet() for k in range(D.poles)],
        "rotor_carrier": [Pos(0, 0, zs + STACK) * rotor_carrier()],
        "jacket_inner": [Pos(0, 0, STATIONS["jacket"]) * jacket_inner()],
        "jacket_outer": [Pos(0, 0, STATIONS["jacket"]) * jacket_outer()],
        "manifold": [Pos(0, 0, STATIONS["manifold"]) * manifold()],
        "drive_plate": [Pos(0, 0, STATIONS["drive_plate"]) * drive_plate()],
        "cap_rear": [Pos(0, 0, STATIONS["cap_rear"]) * cap(front=False)],
        "cap_front": [Pos(0, 0, LENGTH - CAP) * cap(front=True)],
        "bearing": [Pos(0, 0, 0) * bearing(), Pos(0, 0, LENGTH - BRG_W) * bearing()],
        "shaft": [Pos(0, 0, -3.0) * shaft()],
    }
    return parts


def volume_between(a, b):
    try:
        return (a & b).volume
    except Exception:
        return 0.0


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    parts = assembly()
    print(f"motor Ø{D.od:.0f} x {LENGTH:.1f} mm; stack at z = {STATIONS['stack']:.1f}..{STATIONS['stack']+STACK:.1f}; "
          f"end winding overhang {OVERHANG:.1f} mm")
    print(f"coil: {D.turns} turns of {CU_W:.1f} x {CU_H:.1f} mm copper, Ø{D.cond_bore} bore; bend radius "
          f"{bend_radius(0):.1f}..{bend_radius(D.turns-1):.1f} mm on the conductor centre\n")
    print(f"  {'part':14} {'n':>3} {'volume mm³':>11} {'bbox mm':>22}  ok")
    for name, solids in parts.items():
        s = solids[0]
        bb = s.bounding_box().size
        ok = s.is_valid and s.volume > 0
        print(f"  {name:14} {len(solids):3d} {s.volume:11.0f} {bb.X:6.1f} x {bb.Y:5.1f} x {bb.Z:5.1f}  {'ok' if ok else 'INVALID'}")
        export_stl(Compound(children=solids), str(OUT / f"{name}.stl"), tolerance=0.02, angular_tolerance=0.35)
    export_step(Compound(children=[s for v in parts.values() for s in v]), str(OUT / "assembly.step"))

    print("\ninterference (mm³ shared; by symmetry one coil, one magnet):")
    c0, c1 = parts["coil"][0], parts["coil"][1]
    checks = [
        ("coil / stator core", c0, parts["stator_core"][0]),
        ("coil / next coil", c0, c1),
        ("coil / manifold", c0, parts["manifold"][0]),
        ("coil / jacket inner", c0, parts["jacket_inner"][0]),
        ("coil / front cap", c0, parts["cap_front"][0]),
        ("coil / rotor carrier", c0, parts["rotor_carrier"][0]),
        ("magnet / stator core", parts["magnet"][0], parts["stator_core"][0]),
        ("magnet / rotor yoke", parts["magnet"][0], parts["rotor_yoke"][0]),
        ("stator core / jacket inner", parts["stator_core"][0], parts["jacket_inner"][0]),
        ("jacket inner / jacket outer", parts["jacket_inner"][0], parts["jacket_outer"][0]),
        ("rotor carrier / front cap", parts["rotor_carrier"][0], parts["cap_front"][0]),
        ("drive plate / manifold", parts["drive_plate"][0], parts["manifold"][0]),
        ("rotor yoke / manifold", parts["rotor_yoke"][0], parts["manifold"][0]),
    ]
    bad = 0
    for label, a, b in checks:
        v = volume_between(a, b)
        bad += v > 0.01
        print(f"  {label:30} {v:8.3f}  {'ok' if v <= 0.01 else 'CLASH'}")
    # clearances that are design numbers, read back from the solids
    lo = min(leg_xy(0, -1)[0], leg_xy(0, 1)[0])
    print(f"\nconductor to slot wall: {D.liner:.2f} mm liner; turn 0 inner face at r ≈ {turn_radius(0) - CU_H/2:.2f} "
          f"(slot floor {R_SLOT_A:.2f})")
    print(f"{'no clashes' if not bad else str(bad) + ' CLASHES'}; wrote {OUT}/*.stl and assembly.step")
    render()
    return bad


def render():
    """Cutaway PNG through OpenSCAD (cutaway.scad imports the STLs)."""
    import shutil
    import subprocess
    if not (shutil.which("openscad") and shutil.which("xvfb-run")):
        print("render skipped: needs openscad and xvfb-run")
        return
    subprocess.run(["xvfb-run", "-a", "openscad", "-o", "out/cad.png", "--imgsize=1500,1100",
                    "--camera=0,0,34,75,0,95,300", "--colorscheme=Tomorrow", "--projection=p", "cutaway.scad"],
                   cwd=HERE, check=True, capture_output=True)
    print(f"wrote {OUT.parent / 'cad.png'}")


if __name__ == "__main__":
    sys.exit(1 if main() else 0)
