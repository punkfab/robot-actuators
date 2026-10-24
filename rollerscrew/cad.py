"""
CAD of the roller-screw actuator design.py chose, as an assembly, rod retracted.

All dimensions come from design.chosen() (out/design.json) and design's layout
constants; the axial stations are Actuator.stations(), the same numbers the sizing used.
Coaxial, so each part is built in place along z (rear pin at z = 0, rod toward +z).

Load path: rod eye -> screw -> rollers -> nut -> rotor -> front bearing pair -> housing
-> rear cap -> tail tube -> load cell -> rear clevis.

Parts
  screw, rod_eye          threaded rod (the output) and its eye
  roller (xN)             single-start rollers with a sync gear at each end
  ring_gear (x2)          pressed into the nut ends, mesh the roller gears
  carrier (x2)            hold the roller journals at even spacing
  nut                     multi-start internal thread
  rotor, magnet (xP)      hub that carries the nut, yoke over the stack, surface magnets
  stator_core, coil (xS)  concentrated winding, one coil per tooth
  bearing_front (x2), bearing_rear
  housing, cap_front, cap_rear, tail_tube, load_cell, clevis, encoder_magnet, encoder_board, boot

Stand-ins, stated: threads are rings of revolution (no helix), phased so screw, roller
and nut flanks meet; gear teeth are trapezoids, not involutes; coil sides are
rectangular, filling less of the tapered slot than the sizing assumes; bearings are
races and balls with no cage; no fasteners, leads or seals beyond the wiper bore.

    python rollerscrew/cad.py        (or: make rollerscrew-cad)
"""

import sys
from math import cos, degrees, floor, sin, tau
from pathlib import Path

from build123d import (Align, Axis, Box, BuildLine, BuildPart, BuildSketch, Compound, Cylinder, Plane,
                       Polyline, Pos, Rot, Sphere, export_step, export_stl, extrude, make_face, revolve)

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402

OUT = HERE / "out" / "cad"
A = dz.chosen()
S, M, L = A.screw, A.motor, A.link
Z = A.stations()
P = S.pitch
H = P / 4                                   # thread half-height about the pitch radius (45 deg flanks)
R_S, R_R, R_N = S.d_s / 2, S.d_r / 2, S.d_n / 2
R_ORBIT = R_S + R_R
CLEAR = 0.03                                # roller thread shrink, so flanks meet with a gap
PIN_MIN = L.stroke()[1]                     # the linkage's pin-to-pin length, rod retracted
ROD_END = PIN_MIN - dz.EYE

# sync gears: teeth cut into the roller's threaded ends (they keep the thread profile, so
# they still mesh with the screw); the ring gear's teeth reach into the roller's tooth
# spaces. Module about half the pitch; ring teeth = z_r * d_n / d_r.
Z_R = round(S.d_r / (P / 2))
MOD = S.d_r / Z_R
Z_RING = round(Z_R * S.d_n / S.d_r)
GEAR_W = 2.5
JOURNAL_W = dz.NUT_END - GEAR_W

R_HUB_IN = S.nut_od / 2
R_HUB_OUT = R_HUB_IN + 2.0
R_YOKE_OUT = R_HUB_IN + M.rotor_yoke
R_MAG_IN, R_GAP = M.r_gap - M.magnet, M.r_gap
R_SI, R_SO = M.r_si, M.r_so
R_SLOT_A, R_SLOT_B = R_SI + M.tooth_tip, R_SI + M.tooth_tip + M.slot_depth
R_HOUSING = R_SO + 2.0
B_IN = R_HUB_OUT
B_OUT = B_IN + dz.B_SECTION


def ring(r_in, r_out, z0, h):
    c = Cylinder(r_out, h, align=(Align.CENTER, Align.CENTER, Align.MIN))
    if r_in > 0:
        c -= Cylinder(r_in, h, align=(Align.CENTER, Align.CENTER, Align.MIN))
    return Pos(0, 0, z0) * c


def revolved(pts):
    """solid of revolution about z from a closed (r, z) outline"""
    with BuildPart() as p:
        with BuildSketch(Plane.XZ):
            with BuildLine():
                Polyline(*pts, close=True)
            make_face()
        revolve(axis=Axis.Z)
    return p.part


def zigzag(r_p, z0, z1, crest_out, phase, shrink=0.0):
    """thread outline (r, z) from z0 to z1 exactly on a triangle wave of pitch P, 45-degree
    flanks, so a thread starts and ends on its own profile. crest_out: crests point away
    from the axis (external thread). phase: z of a crest, mod P."""
    def r(z):
        t = ((z - phase) / P) % 1.0
        h = H - 4 * H * min(t, 1 - t)                      # +H at a crest, -H half a pitch away
        return (r_p + h if crest_out else r_p - h) - shrink
    zs = [z0] + [phase + k * P / 2 for k in range(floor((z0 - phase) / (P / 2)) + 1, floor((z1 - phase) / (P / 2)) + 1)
                 if z0 < phase + k * P / 2 < z1] + [z1]
    return [(r(z), z) for z in zs]


def gear_outline(z_teeth, tip, root):
    """trapezoid teeth between a tip and a root radius (tip < root for a ring gear)"""
    pts, w = [], tau / z_teeth
    for i in range(z_teeth):
        a = i * w
        for da, r in ((-0.30, root), (-0.16, tip), (0.16, tip), (0.30, root)):
            pts.append((r * cos(a + da * w), r * sin(a + da * w)))
    return pts


def prism(pts, z0, h):
    with BuildPart() as p:
        with BuildSketch(Plane.XY.offset(z0)):
            with BuildLine():
                Polyline(*pts, close=True)
            make_face()
        extrude(amount=h)
    return p.part


# ---- the roller screw -----------------------------------------------------------------
Z_NT0 = Z["nut_thread"]
Z_NT1 = Z_NT0 + S.nut_len
Z_NUT0, Z_NUT1 = Z["nut"], Z["nut"] + S.nut_len + 2 * dz.NUT_END


def screw():
    z0, z1 = Z["screw_rear"], ROD_END
    pts = [(0, z0)] + zigzag(R_S, z0, z1, True, 0.0) + [(0, z1)]
    return revolved(pts)


def rod_eye():
    w, h = 10.0, dz.EYE + 8.0
    eye = Pos(0, 0, ROD_END) * Box(w, S.d_s + 4, h, align=(Align.CENTER, Align.CENTER, Align.MIN))
    eye += Pos(0, 0, PIN_MIN) * Rot(0, 90, 0) * Cylinder((S.d_s + 4) / 2, w)
    return eye - Pos(0, 0, PIN_MIN) * Rot(0, 90, 0) * Cylinder(3.0, w + 2)


ROLLER_ROOT = R_R - H - 0.12               # roller tooth-space bottom, below the thread root


def roller_local():
    """one roller on its own axis (z global): thread over the nut's threaded length and on
    through both gear zones, tooth spaces cut into the gear zones, plain journals beyond"""
    r_j = R_R * 0.5
    zg0, zg1 = Z_NT0 - GEAR_W, Z_NT1 + GEAR_W
    zj0, zj1 = Z_NT0 - dz.NUT_END, Z_NT1 + dz.NUT_END
    pts = [(0, zj0), (r_j, zj0), (r_j, zg0)] + zigzag(R_R, zg0, zg1, True, P / 2, CLEAR) + [(r_j, zg1), (r_j, zj1), (0, zj1)]
    body = revolved(pts)
    star = gear_outline(Z_R, R_R + H + 0.5, ROLLER_ROOT)
    for z0 in (zg0, Z_NT1):
        spaces = ring(ROLLER_ROOT, R_R + H + 0.4, z0, GEAR_W) - prism(star, z0 - 0.1, GEAR_W + 0.2)
        body -= spaces
    return body


def roller_placements():
    """(orbit angle, own rotation) for each roller: own rotation puts a tooth gap where
    the ring has a tooth tip at that roller's mesh point"""
    out = []
    pitch = tau * R_N / Z_RING
    for k in range(S.n_rollers):
        th = tau * k / S.n_rollers
        j = round(th * Z_RING / tau)
        s_tip = R_N * (tau * j / Z_RING - th)                # ring tip, arc distance from the mesh point
        beta = th + (s_tip + pitch / 2) / R_R                # roller tooth tip half a pitch from the ring's
        out.append((th, beta))
    return out


def rollers():
    base = roller_local()
    return [Pos(R_ORBIT * cos(th), R_ORBIT * sin(th), 0) * Rot(0, 0, degrees(beta)) * base
            for th, beta in roller_placements()]


def ring_gears():
    """tips reach into the roller tooth spaces (inside the roller thread band), roots clear the roller crests"""
    out = []
    star = gear_outline(Z_RING, R_N - H + 0.04, R_N + H + 0.15)
    for z0 in (Z_NT0 - GEAR_W, Z_NT1):
        outer = ring(0, S.nut_od / 2 - 2.5, z0, GEAR_W)            # sits in the nut's end counterbore
        out.append(outer - prism(star, z0 - 0.1, GEAR_W + 0.2))
    return out


def carriers():
    r_in, r_out = R_S + H + 0.4, R_N + H + 0.1
    out = []
    for z0 in (Z_NUT0, Z_NT1 + GEAR_W):
        c = ring(r_in, r_out, z0, JOURNAL_W)
        for th, _ in roller_placements():
            c -= Pos(R_ORBIT * cos(th), R_ORBIT * sin(th), z0 - 0.1) * Cylinder(R_R * 0.5 + 0.05, JOURNAL_W + 0.2,
                                                                                 align=(Align.CENTER, Align.CENTER, Align.MIN))
        out.append(c)
    return out


def nut():
    """internal thread over the threaded length; counterbores at each end take a ring
    gear and a carrier"""
    r_o, cb = S.nut_od / 2, S.nut_od / 2 - 2.5
    pts = [(r_o, Z_NUT0), (cb, Z_NUT0), (cb, Z_NT0)] + zigzag(R_N, Z_NT0, Z_NT1, False, 0.0) + [(cb, Z_NT1), (cb, Z_NUT1), (r_o, Z_NUT1)]
    return revolved(pts)


# ---- motor ----------------------------------------------------------------------------
def rotor():
    hub = ring(R_HUB_IN, R_HUB_OUT, Z["encoder"], Z["cap_front"] - Z["encoder"])
    return hub + ring(R_HUB_IN, R_YOKE_OUT, Z["stack"], M.stack)


def magnets():
    arc = 0.85 * tau / M.poles
    out = []
    for i in range(M.poles):
        a0 = i * tau / M.poles - arc / 2
        n = 8
        pts = [(R_MAG_IN * cos(a0 + arc * t / n), R_MAG_IN * sin(a0 + arc * t / n)) for t in range(n + 1)]
        pts += [(R_GAP * cos(a0 + arc * t / n), R_GAP * sin(a0 + arc * t / n)) for t in range(n, -1, -1)]
        out.append(prism(pts, Z["stack"], M.stack))
    return out


def stator_core():
    core = ring(R_SI, R_SO, Z["stack"], M.stack) - ring(R_SLOT_A, R_SLOT_B, Z["stack"] - 1, M.stack + 2)
    for i in range(M.slots):
        a = degrees(i * tau / M.slots)
        core += Rot(0, 0, a) * Pos((R_SI + R_SLOT_B) / 2 + 0.05, 0, Z["stack"]) * Box(
            R_SLOT_B - R_SI + 0.1, M.tooth_w, M.stack, align=(Align.CENTER, Align.CENTER, Align.MIN))
    open_w = 0.2 * tau * R_SI / M.slots
    for i in range(M.slots):
        a = degrees((i + 0.5) * tau / M.slots)
        core -= Rot(0, 0, a) * Pos((R_SI + R_SLOT_A) / 2, 0, Z["stack"] - 1) * Box(
            M.tooth_tip + 1.0, open_w, M.stack + 2, align=(Align.CENTER, Align.CENTER, Align.MIN))
    return core


COIL_T = (tau * R_SLOT_A / M.slots - M.tooth_w) / 2 - 0.35          # coil side width, set by the slot's narrow end


def coil_local():
    """one tooth's coil, tooth on +x"""
    y_max = M.tooth_w / 2 + 0.2 + COIL_T
    x0, x1 = R_SLOT_A + 0.2, (( R_SLOT_B - 0.2) ** 2 - y_max ** 2) ** 0.5
    rl, xc = x1 - x0, (x0 + x1) / 2
    zl = M.stack + 2 * (COIL_T + 0.3)
    outer = Pos(xc, 0, Z["stack"] - COIL_T - 0.3) * Box(rl, M.tooth_w + 0.4 + 2 * COIL_T, zl, align=(Align.CENTER, Align.CENTER, Align.MIN))
    inner = Pos(xc, 0, Z["stack"] - 0.3) * Box(rl + 2, M.tooth_w + 0.4, M.stack + 0.6, align=(Align.CENTER, Align.CENTER, Align.MIN))
    return outer - inner


def coils():
    c = coil_local()
    return [Rot(0, 0, degrees(i * tau / M.slots)) * c for i in range(M.slots)]


# ---- bearings, housing, ends -------------------------------------------------------------
def bearing(z0):
    sec = dz.B_SECTION
    inner = ring(B_IN, B_IN + sec * 0.3, z0, dz.B_W)
    outer = ring(B_OUT - sec * 0.3, B_OUT, z0, dz.B_W)
    rb = (B_IN + B_OUT) / 2
    d_ball = sec * 0.42
    n = floor(tau * rb / (d_ball * 1.25))
    balls = [Pos(rb * cos(i * tau / n), rb * sin(i * tau / n), z0 + dz.B_W / 2) * Sphere(d_ball / 2) for i in range(n)]
    return [inner, outer] + balls


def housing():
    z0, z1 = Z["cap_rear"] + dz.CAP, Z["cap_front"]
    h = ring(R_SO, R_HOUSING, z0, z1 - z0)
    h += ring(B_OUT, R_SO, Z["bearing_rear"], dz.B_W)                 # rear bearing seat web
    h += ring(B_OUT, R_SO, Z["bearing_front"], 2 * dz.B_W)            # front bearing pair seat
    return h


def cap_front():
    return ring(R_S + H + 0.2, R_HOUSING, Z["cap_front"], dz.CAP)     # wiper bore


def cap_rear():
    return ring(R_S + H + 0.5, R_HOUSING, Z["cap_rear"], dz.CAP)


def tail_tube():
    r_in = R_S + H + 0.5
    return ring(r_in, r_in + 2.0, Z["tail"], Z["tail_len"] + 0.0) if Z["tail_len"] > 0 else None


def load_cell():
    return ring(0, 12.0, Z["load_cell"], dz.LOAD_CELL)


def clevis():
    w = 10.0
    c = Pos(0, 0, -8) * Box(w * 2.4, 16, dz.CLEVIS + 8, align=(Align.CENTER, Align.CENTER, Align.MIN))
    c -= Pos(0, 0, -9) * Box(w + 0.6, 18, 13, align=(Align.CENTER, Align.CENTER, Align.MIN))     # fork slot
    return c - Rot(0, 90, 0) * Cylinder(3.0, w * 3)


def encoder():
    mag = ring(R_HUB_OUT, R_HUB_OUT + 1.5, Z["encoder"] + 0.5, 2.0)
    board = ring(R_HUB_OUT + 2.5, R_SO - 0.2, Z["encoder"] + 3.5, 1.6)
    return mag, board


def boot():
    """bellows over the rod between the front cap and the eye, compressed (rod retracted)"""
    z0, z1 = Z["front"], ROD_END
    r_a, r_b = R_S + 1.5, R_S + 4.0
    n = max(3, int((z1 - z0) / 2.5))
    pts = [(r_a - 0.6, z0)]
    for i in range(n + 1):
        pts.append((r_a if i % 2 == 0 else r_b, z0 + (z1 - z0) * i / n))
    pts += [(r_a - 0.6, z1)]
    inner = [(r - 0.6, z) for r, z in pts[1:-1]][::-1]
    return revolved(pts[:-1] + [pts[-1]] + inner[:-1]) if False else ring(R_S + H + 0.3, r_a, z0, z1 - z0) + revolved(
        [(r_a - 0.01, z0)] + pts[1:-1] + [(r_a - 0.01, z1)])


def assembly():
    parts = {
        "screw": [screw()], "rod_eye": [rod_eye()], "roller": rollers(), "ring_gear": ring_gears(),
        "carrier": carriers(), "nut": [nut()], "rotor": [rotor()], "magnet": magnets(),
        "stator_core": [stator_core()], "coil": coils(),
        "bearing_front": bearing(Z["bearing_front"]) + bearing(Z["bearing_front"] + dz.B_W),
        "bearing_rear": bearing(Z["bearing_rear"]),
        "housing": [housing()], "cap_front": [cap_front()], "cap_rear": [cap_rear()],
        "load_cell": [load_cell()], "clevis": [clevis()], "encoder_magnet": [encoder()[0]],
        "encoder_board": [encoder()[1]], "boot": [boot()],
    }
    t = tail_tube()
    if t is not None:
        parts["tail_tube"] = [t]
    return parts


def shared(a, b):
    ba, bb = a.bounding_box(), b.bounding_box()
    if any(tuple(ba.min)[i] > tuple(bb.max)[i] or tuple(bb.min)[i] > tuple(ba.max)[i] for i in range(3)):
        return 0.0
    try:
        r = a.intersect(b)
    except Exception as e:                      # OCC sometimes fails on coincident faces
        print(f"    (intersection failed: {type(e).__name__})")
        return float("nan")
    return r.volume if r is not None else 0.0


CHECKS = [("screw", "roller"), ("roller", "nut"), ("roller", "roller"), ("roller", "ring_gear"), ("roller", "carrier"),
          ("screw", "carrier"), ("screw", "nut"), ("screw", "ring_gear"), ("nut", "rotor"), ("ring_gear", "nut"),
          ("magnet", "stator_core"), ("rotor", "stator_core"), ("coil", "stator_core"), ("coil", "coil"),
          ("coil", "rotor"), ("coil", "housing"), ("bearing_front", "housing"), ("bearing_front", "rotor"),
          ("bearing_rear", "rotor"), ("stator_core", "housing"), ("screw", "cap_front"), ("screw", "cap_rear"),
          ("screw", "rotor"), ("encoder_board", "rotor"), ("encoder_magnet", "encoder_board"), ("boot", "screw")]


def interference(parts):
    print("\ninterference, mm^3 shared (0 = clear or face contact):")
    worst = {}
    for a, b in CHECKS:
        va = parts[a] if a != "roller" else parts[a][:2]
        vb = parts[b] if b != "roller" else parts[b][:2]
        if a == b:
            pairs = [(va[0], va[1])] if len(va) > 1 else []
        else:
            pairs = [(x, y) for x in va for y in vb]
        v = sum(shared(x, y) for x, y in pairs)
        worst[f"{a}/{b}"] = v
        print(f"  {a:>15} / {b:<15} {v:10.4f}")
    return worst


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"screw Ø{S.d_s} x {S.lead} lead, {S.starts} starts, pitch {P}; {S.n_rollers} rollers Ø{S.d_r:.2f}; "
          f"nut thread {S.nut_len:.1f} mm; sync gears m{MOD:.2f}, {Z_R}/{Z_RING} teeth")
    print(f"motor Ø{M.od:.1f} x {M.stack} stack; pins {PIN_MIN:.1f} mm retracted, stroke {Z['stroke']:.1f} mm")
    parts = assembly()
    for name, solids in parts.items():
        export_stl(Compound(children=solids), str(OUT / f"{name}.stl"), tolerance=0.01, angular_tolerance=0.2)
        bb = Compound(children=solids).bounding_box()
        print(f"  {name:15} x{len(solids):<3} z {bb.min.Z:7.1f} .. {bb.max.Z:7.1f}  r <= {max(abs(bb.max.X), abs(bb.min.X), abs(bb.max.Y)):5.1f}")
    export_step(Compound(children=[s for v in parts.values() for s in v]), str(OUT / "assembly.step"))
    interference(parts)
    print(f"wrote {OUT}/*.stl and assembly.step; renders: make rollerscrew-render")


if __name__ == "__main__":
    main()
