"""
Parametric CAD for the 2-CELL VARIABLE CHAIN MOTOR (see chain.py for the physics).

Two Goolsky 2204 outrunners (Ø28 x 22 mm, 3 mm shaft exiting the MOUNT side — see the root
README), each turned into a motor cell by a STATOR GEAR (under the mount face; its hollow hub
clears the shaft) and a ROTOR CAP (gear + skirt bonded over the bell + a 5 mm steel dowel on top
for the upper bearing). Identical m1 z34 spur gears on both, so r_r = r_s = 1 and the link's
centre distance is one pitch diameter (34 mm).

    base          fixed to the "upper arm"; cell A's hub keys into a hex pocket + one M3
    link_lower    swing link, lower plate: 6800 (10x19x5) seats around both stator hubs —
                  pivots on A's hub, carries B's hub
    link_upper    swing link, upper plate + 2 fused posts: 685 (5x11x5) seats on both rotor-cap
                  dowels; A's dowel continues up as the OUTPUT (encoder magnet on top)
    stator_gear_A gear + LONG hub (through link_lower, hex end into the base) — never turns
    stator_gear_B gear + SHORT hub (spins in link_lower) — turns 2φ in the world as B orbits
    rotor_gear    x2, "rotor cap": gear + skirt bonded over the bell + press hole for the dowel
    dowel         reference only: Ø5 steel dowel pin (rotor cap -> 685 bearing / output)
    motor_2204    reference body only (not printed): mount + bell + 3 mm shaft out the bottom

Z stack (mm):  base [-9,-1] | link_lower [0,5] | stator gears [5.5,9.5] | motor [9.5,31.5]
               | rotor caps: skirt [25.5,31.5] + gear [31.5,35.5] | link_upper [37,42]

Interference: `check_interference()` intersects every pair of placed solids at several hinge
angles (idea borrowed from earthtojake/text-to-cad's overlap_volume checks) — the gear meshes
must phase correctly (B offset by half a tooth) for it to pass.

    ../.venv/bin/python chain-motor/cad.py       # parts -> out/*.step/*.stl + interference report
"""

import sys
from dataclasses import dataclass
from math import pi, cos, sin, degrees
from pathlib import Path

from build123d import (
    Box, Cylinder, Pos, Rot, Compound, Align, RegularPolygon, extrude,
    export_step, export_stl,
)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "cycloidal"))
from chain import ChainParams
from gears import spur_gear

_BOT = (Align.CENTER, Align.CENTER, Align.MIN)


@dataclass
class Cad:
    P: ChainParams

    # --- 2204 envelope (reference) ---
    bell_d: float = 28.0
    motor_h: float = 22.0
    shaft_d: float = 3.0              # exits the mount side (down, through the stator hub)
    dowel_d: float = 5.0              # steel dowel in the rotor cap -> 685 bearing / output
    skirt_h: float = 6.0              # rotor-cap skirt over the bell top
    skirt_wall: float = 1.5
    mount_pat: tuple = (16.0, 19.0)   # M3 cross pattern on the stator mount
    m3: float = 3.3

    # --- stack ---
    gear_t: float = 4.0
    plate_t: float = 5.0
    gap: float = 0.5
    b6800: tuple = (10.0, 19.0, 5.0)  # hub bearing ID, OD, W
    b685: tuple = (5.0, 11.0, 5.0)    # shaft bearing
    link_r: float = 20.0              # link plate end radius
    post_r: float = 3.0
    post_y: float = 15.0              # posts at (a/2, ±post_y), clear of both gear tips
    base_r: float = 22.0
    hex_r: float = 4.95               # A-hub anti-rotation hex (circumradius; fits the Ø10 journal)
    m3_tap: float = 2.5

    @property
    def a(self): return self.P.center_mm
    @property
    def z_stator(self): return self.plate_t + self.gap                      # 5.5
    @property
    def z_motor(self): return self.z_stator + self.gear_t                   # 9.5
    @property
    def z_rotor(self): return self.z_motor + self.motor_h                   # 28.5
    @property
    def z_upper(self): return self.z_rotor + self.gear_t + 1.5              # 34
    @property
    def tooth_phase_deg(self): return 180.0 / self.P.z_rotor                 # half a tooth


def _mount_holes(c: Cad, h: float, z0: float):
    cut = None
    for i, r in enumerate(c.mount_pat):
        for k in (0, 1):
            ang = pi / 2 * (2 * k + i) + pi / 4           # opposite pairs on the two diagonals
            hole = Pos(r / 2 * cos(ang), r / 2 * sin(ang), z0) * Cylinder(c.m3 / 2, h, align=_BOT)
            cut = hole if cut is None else cut + hole
    return cut


def make_stator_gear(c: Cad, long_hub: bool):
    """Stator gear: m1 z34 + hollow hub down into the link's 6800 bearing.

    A (long_hub): below the link the hub turns into a HEX that keys into the base (it carries
    the whole stator reaction torque), with a central M3 tapped hole for axial retention. The
    Ø6 shaft/circlip bore stops above the tap. (First cut used 3x M2 at r=3.9 in a Ø6-bore,
    Ø10 hub — the holes broke into the bore: 0.11 mm wall, flagged by text-to-cad's dfam-check.)
    """
    P = c.P
    g = Pos(0, 0, c.z_stator) * spur_gear(P.module_mm, P.z_stator, c.gear_t)
    hub_id = 6.0                                     # clears the 2204's shaft circlip
    if long_hub:
        hub = Pos(0, 0, -1.0) * Cylinder(c.b6800[0] / 2, c.z_stator + 1.01, align=_BOT)
        hub += Pos(0, 0, -7.0) * extrude(RegularPolygon(c.hex_r, 6), amount=6.01)
        part = g + hub
        part -= Pos(0, 0, 2.0) * Cylinder(hub_id / 2, 30, align=_BOT)
        part -= Pos(0, 0, -8.0) * Cylinder(c.m3_tap / 2, 10.01, align=_BOT)
    else:
        hub = Cylinder(c.b6800[0] / 2, c.z_stator + 0.01, align=_BOT)
        part = g + hub
        part -= Pos(0, 0, -1) * Cylinder(hub_id / 2, 30, align=_BOT)
    part -= _mount_holes(c, c.gear_t + 2, c.z_stator - 1)
    return part


def make_rotor_gear(c: Cad):
    """Rotor cap: m1 z34 gear on top of a skirt bonded over the bell, with a press-fit hole for
    a Ø5 steel dowel (up through the link's 685 bearing). Prints gear-face-down, support-free.
    (First cut fused the dowel into the print — 13-22 % support in every orientation.)"""
    P = c.P
    g = Pos(0, 0, c.z_rotor) * spur_gear(P.module_mm, P.z_rotor, c.gear_t)
    r_in = c.bell_d / 2 + 0.1
    skirt = Pos(0, 0, c.z_rotor - c.skirt_h) * Cylinder(r_in + c.skirt_wall, c.skirt_h + 0.01, align=_BOT)
    skirt -= Pos(0, 0, c.z_rotor - c.skirt_h - 1) * Cylinder(r_in, c.skirt_h + 1, align=_BOT)
    part = g + skirt
    part -= Pos(0, 0, c.z_rotor - 1) * Cylinder(c.dowel_d / 2 - 0.05, c.gear_t + 2, align=_BOT)
    return part


def make_dowel(c: Cad):
    """Reference only (not printed): Ø5 steel dowel pressed into the rotor cap."""
    top = c.z_upper + c.plate_t + 3.0
    return Pos(0, 0, c.z_rotor) * Cylinder(c.dowel_d / 2 - 0.06, top - c.z_rotor, align=_BOT)


def _stadium(c: Cad, r: float, h: float, z0: float):
    body = Pos(0, 0, z0) * Cylinder(r, h, align=_BOT)
    body += Pos(c.a, 0, z0) * Cylinder(r, h, align=_BOT)
    body += Pos(c.a / 2, 0, z0) * Box(c.a, 2 * r, h, align=_BOT)
    return body


def make_link_lower(c: Cad):
    p = _stadium(c, c.link_r, c.plate_t, 0.0)
    for x in (0.0, c.a):                             # 6800 seat + hub through-bore
        p -= Pos(x, 0, -1) * Cylinder(c.b6800[1] / 2 + 0.1, c.plate_t + 2, align=_BOT)
    for s in (1, -1):                                # M3 tapped holes for the posts' screws
        p -= Pos(c.a / 2, s * c.post_y, -1) * Cylinder(1.25, c.plate_t + 2, align=_BOT)
    return p


def make_link_upper(c: Cad):
    t = c.plate_t
    p = _stadium(c, 9.0, t, c.z_upper)               # slim: only the shaft bearings up here
    p += Pos(c.a / 2, 0, c.z_upper) * Box(8.0, 2 * c.post_y, t, align=_BOT)
    for s in (1, -1):                                # posts down to link_lower
        p += Pos(c.a / 2, s * c.post_y, c.plate_t) * Cylinder(c.post_r, c.z_upper - c.plate_t + 0.01,
                                                              align=_BOT)
        p -= Pos(c.a / 2, s * c.post_y, c.plate_t - 1) * Cylinder(c.m3 / 2 - 0.5, 12, align=_BOT)
    for x in (0.0, c.a):                             # 685 seats
        p -= Pos(x, 0, c.z_upper - 1) * Cylinder(c.b685[1] / 2 + 0.1, t + 2, align=_BOT)
    return p


def make_base(c: Cad):
    """Round flange under A (+ a tab toward -x to bolt to the upper arm). A hex pocket takes
    A's hub (torque), an M3 through the floor pulls it down (retention)."""
    p = Pos(0, 0, -9.0) * Cylinder(c.base_r * 0.55, 8.0, align=_BOT)
    p += Pos(-c.base_r, 0, -9.0) * Box(2 * c.base_r, 16.0, 8.0, align=_BOT)
    p -= Pos(0, 0, -7.0) * extrude(RegularPolygon(c.hex_r + 0.15, 6), amount=7)
    p -= Pos(0, 0, -10) * Cylinder(c.m3 / 2, 4, align=_BOT)
    for x in (-c.base_r - 8, -c.base_r + 8):
        p -= Pos(x - 4, 0, -10) * Cylinder(c.m3 / 2, 10, align=_BOT)
    return p


def make_motor(c: Cad):
    """Reference 2204 envelope: stator mount + bell, 3 mm shaft out the MOUNT side (down)."""
    mount = Pos(0, 0, c.z_motor) * Cylinder(12.0, 4.0, align=_BOT)
    bell = Pos(0, 0, c.z_motor + 4.0) * Cylinder(c.bell_d / 2, c.motor_h - 4.0, align=_BOT)
    shaft = Pos(0, 0, c.z_motor - 6.0) * Cylinder(c.shaft_d / 2, 6.01, align=_BOT)
    return mount + bell + shaft


def build_all(c: Cad) -> dict:
    return {
        "base": make_base(c),
        "link_lower": make_link_lower(c),
        "link_upper": make_link_upper(c),
        "stator_gear_A": make_stator_gear(c, True),
        "stator_gear_B": make_stator_gear(c, False),
        "rotor_gear": make_rotor_gear(c),
        "dowel": make_dowel(c),
        "motor_2204": make_motor(c),
    }


def placed(c: Cad, parts: dict, phi_deg=0.0, theta_deg=0.0) -> dict:
    """World placement at hinge φ and rotor θ (chain.py kinematics, r_r = r_s = 1):
       A stator 0, A rotor θ; link φ; B stator 2φ (planet on a fixed sun), B rotor 2φ − θ;
       B's gears carry the extra half-tooth phase so the teeth interleave."""
    ph = c.tooth_phase_deg
    B = Rot(0, 0, phi_deg) * Pos(c.a, 0, 0)       # B's local frame is already turned by φ,
                                                  # so a world angle of 2φ is a local φ
    return {
        "base": parts["base"],
        "link_lower": Rot(0, 0, phi_deg) * parts["link_lower"],
        "link_upper": Rot(0, 0, phi_deg) * parts["link_upper"],
        "stator_gear_A": parts["stator_gear_A"],
        "motor_A": parts["motor_2204"],
        "rotor_gear_A": Rot(0, 0, theta_deg) * parts["rotor_gear"],
        "dowel_A": parts["dowel"],
        "stator_gear_B": B * Rot(0, 0, phi_deg + ph) * parts["stator_gear_B"],
        "motor_B": B * Rot(0, 0, phi_deg) * parts["motor_2204"],
        "rotor_gear_B": B * Rot(0, 0, phi_deg - theta_deg + ph) * parts["rotor_gear"],
        "dowel_B": B * parts["dowel"],
    }


# intended contacts that the interference check should ignore (shaft/hub through bores etc.)
_ALLOWED = {
    frozenset(("stator_gear_A", "base")), frozenset(("stator_gear_A", "link_lower")),
    frozenset(("stator_gear_B", "link_lower")), frozenset(("link_upper", "link_lower")),
    frozenset(("motor_A", "stator_gear_A")), frozenset(("motor_B", "stator_gear_B")),
    frozenset(("motor_A", "rotor_gear_A")), frozenset(("motor_B", "rotor_gear_B")),
    frozenset(("dowel_A", "rotor_gear_A")), frozenset(("dowel_B", "rotor_gear_B")),
}
_MESHES = (("rotor_gear_A", "rotor_gear_B"), ("stator_gear_A", "stator_gear_B"))


def check_interference(c: Cad, parts: dict, phis=(0, 45, 90, 155, -120), thetas=(0, 7.3)):
    """Pairwise intersection volume of every placed solid, plus the minimum CLEARANCE across each
    gear mesh (overlap volume alone reports touching teeth as 0 — zero backlash slips through).
    Returns (worst unexpected overlap mm³, min gear-mesh clearance mm)."""
    worst, gap = 0.0, float("inf")
    for phi in phis:
        for th in thetas:
            pl = placed(c, parts, phi, th)
            names = list(pl)
            for i in range(len(names)):
                for j in range(i + 1, len(names)):
                    a, b = names[i], names[j]
                    if frozenset((a, b)) in _ALLOWED:
                        continue
                    bba, bbb = pl[a].bounding_box(), pl[b].bounding_box()
                    if (bba.min.X > bbb.max.X or bbb.min.X > bba.max.X or bba.min.Y > bbb.max.Y
                            or bbb.min.Y > bba.max.Y or bba.min.Z > bbb.max.Z or bbb.min.Z > bba.max.Z):
                        continue
                    x = pl[a] & pl[b]
                    v = 0.0 if x is None else sum(s.volume for s in (x.solids() if hasattr(x, "solids") else x))
                    if v > 1e-3:
                        print(f"   ✗ φ={phi:4}° θ={th:4}°  {a} ∩ {b} = {v:.3f} mm³")
                    worst = max(worst, v)
            if phi in phis[:2] and th == thetas[0]:     # clearance is φ-invariant; it's slow
                for a, b in _MESHES:
                    gap = min(gap, pl[a].distance_to(pl[b]))
    return worst, gap


def validate(c: Cad):
    tip = c.P.module_mm * (c.P.z_rotor / 2 + 1)
    checks = {
        "bells clear each other (a > bell Ø)": c.a > c.bell_d + 1.0,
        "posts clear both gear tips": ((c.a / 2) ** 2 + c.post_y ** 2) ** 0.5 > tip + c.post_r + 0.5,
        "stator hub fits 6800 bore": abs(c.b6800[0] - 10.0) < 1e-6,
        "dowel fits 685 bore": abs(c.b685[0] - c.dowel_d) < 1e-6,
        "A-hub hex fits the Ø10 journal": c.hex_r <= c.b6800[0] / 2,
        "A-hub wall around M3 tap ≥ 1.2 mm (FDM)": c.hex_r * cos(pi / 6) - c.m3_tap / 2 >= 1.2,
        "rotor-cap skirt wall ≥ 1.2 mm (FDM)": c.skirt_wall >= 1.2,
        "rotor-cap skirts clear each other": 2 * (c.bell_d / 2 + 0.1 + c.skirt_wall) < c.a - 1.0,
        "stator hub clears the 3 mm shaft": 6.0 > c.shaft_d + 1.5,
        "r_r == r_s (hinge commutation-neutral)": c.P.z_rotor == c.P.z_stator,
    }
    for k, v in checks.items():
        print(f"   {'✓' if v else '✗'} {k}")
    return all(checks.values())


if __name__ == "__main__":
    import trimesh
    c = Cad(ChainParams())
    outdir = HERE / "out"
    outdir.mkdir(exist_ok=True)
    parts = build_all(c)
    print("parts:")
    for name, part in parts.items():
        export_step(part, str(outdir / f"{name}.step"))
        export_stl(part, str(outdir / f"{name}.stl"))
        m = trimesh.load(outdir / f"{name}.stl")
        bodies = len(m.split(only_watertight=False))
        print(f"   {name:14s} bodies={bodies}  watertight={m.is_watertight}  "
              f"bbox={(m.bounds[1]-m.bounds[0]).round(1)}")
    print("design rules:")
    ok = validate(c)
    print("interference (all pairs, φ ∈ {0,45,90,155,-120}°, θ ∈ {0,7.3}°):")
    worst, gap = check_interference(c, parts)
    print(f"   worst unexpected overlap {worst:.4f} mm³  {'✓' if worst <= 1e-3 else '✗'}")
    print(f"   min gear-mesh clearance {gap:.3f} mm  {'✓' if gap > 0.02 else '✗ (teeth touching: no backlash)'}")
    ok = ok and gap > 0.02
    export_step(Compound(children=list(placed(c, parts, 60, 0).values())), str(outdir / "assembly.step"))
    print(f"wrote {outdir}/*.step, *.stl, assembly.step (φ=60°)  — {'VALID' if ok and worst <= 1e-3 else 'CHECK'}")
