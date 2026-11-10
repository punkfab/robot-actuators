"""
CAD of wobble2 (design.py): one wobble ring between an input carrier with three ball bearings and
an output carrier with twelve stud-type cam followers, placed by the reducer's own kinematics.

    motor | housing_in (spherical seat, lower half) | ring | housing_out (upper half) | output flange
           input carrier: 3 x 6000 on studs at -24 deg     output carrier: 12 x CF4R at +24 deg

The ring is a band cut from a sphere: its outer surface sits in the housing's spherical bore (split
at the equator to assemble), two pins at its equator run in slots that let it nutate but not spin,
and its two end faces are the waves. The faces come from design.py's roller paths, offset at each
radius along the bearing by that radius's own roller angle asin(D / 2 rho) (the rollers are
cylinders, not cones, so the face is not quite ruled through the centre). They are built as a mesh
(manifold3d booleans, the build123d-part fallback): ~2900 azimuth stations hold the wave. The pin
slots are the swept volume of the pins over a nutation cycle. Everything else is B-rep.

    python wobble2/cad.py      (or: make wobble2-cad)   -> out/cad/*.stl, parts.json, checks
"""

import json
import sys
from pathlib import Path

import numpy as np
import trimesh
from build123d import Axis, Box, Cylinder, Polyline, Pos, export_stl, make_face, revolve

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "cad"
sys.path.insert(0, str(HERE))
import design as D  # noqa: E402

W = D.Wobble2()

# ---------------------------------------------------------------- dimensions (mm)
RHO_IN = W.R - 11.0                      # ring inner sphere
RHO_SEAT = D.seat_radius(W)              # ring outer sphere = housing bore
SEAT_GAP = 0.05
WALL = 8.0                               # housing wall beyond the cavity
R_OUT = RHO_SEAT + 9.0                   # housing outer radius
PIN_D, PIN_IN, PIN_OUT = 4.0, 4.0, 5.0   # anti-rotation pins: diameter, depth in the ring, reach into the housing
SHAFT_D = 12.0
SB = dict(name="6001-2RS", d=12.0, D=28.0, B=8.0)            # input shaft bearing in the housing floor
OB = dict(name="61808-2RS (pair)", d=40.0, D=52.0, B=7.0)    # output carrier main bearings
MOTOR = dict(name="pancake BLDC (stand-in)", D=84.0, L=28.0)
WEB_T = 7.0                              # carrier web thickness


def sph(rho, el, az):
    return np.stack([rho * np.cos(el) * np.cos(az), rho * np.cos(el) * np.sin(az), rho * np.sin(el)], -1)


def to_matrix(Rm, o=(0, 0, 0)):
    M = np.eye(4)
    M[:3, :3] = Rm
    M[:3, 3] = o
    return M


def rot_to(u):
    z = np.array([0, 0, 1.0])
    ax = np.cross(z, u)
    s = np.linalg.norm(ax)
    if s < 1e-12:
        return np.eye(3)
    return trimesh.transformations.rotation_matrix(np.arctan2(s, z @ u), ax / s)[:3, :3]


# ---------------------------------------------------------------- the ring (mesh)

def face_curves(row, rhos, n_az=2880, samples=9000):
    """elevation (rad) of this row's face on a uniform azimuth grid, ring frame, at each radius
    rho: the roller-axis path offset by asin(D / 2 rho) (clamped to the bearing's own width)"""
    u = D.axes(row)[:1]
    span = 2 * np.pi * (abs(W.ratio) if row is W.out else 1)
    phis = np.linspace(0, span, samples, endpoint=False)
    P, N = [], []
    for q in phis:
        p, t, n, k, sp = D.frame(W, row, q, u)
        P.append(p[0])
        N.append(n[0])
    P, N = np.array(P), np.array(N)
    az_grid = np.linspace(0, 2 * np.pi, n_az, endpoint=False)
    out = []
    for rho in rhos:
        rc = np.clip(rho, W.R - row.bearing.B / 2, W.R + row.bearing.B / 2)
        g = np.arcsin(row.bearing.D / 2 / rc) * (rc / rho if rho > rc else 1.0)
        pts = np.cos(g) * P - row.side * np.sin(g) * N
        a = np.mod(np.arctan2(pts[:, 1], pts[:, 0]), 2 * np.pi)
        e = np.arcsin(np.clip(pts[:, 2], -1, 1))
        o = np.argsort(a)
        a, e = a[o], e[o]
        a = np.concatenate([a[-1:] - 2 * np.pi, a, a[:1] + 2 * np.pi])
        e = np.concatenate([e[-1:], e, e[:1]])
        out.append(np.interp(az_grid, a, e))
    return az_grid, np.array(out)


def ring_mesh():
    """the band between the front face (input rollers below it) and the back face (output rollers
    above it), RHO_IN..RHO_SEAT, in the ring frame"""
    rhos = np.unique(np.concatenate([[RHO_IN], np.linspace(W.R - 4, W.R + 4, 9), [RHO_SEAT]]))
    az, ef = face_curves(W.inp, rhos)
    _, eb = face_curves(W.out, rhos)
    n, m = len(az), len(rhos)
    # section polygon per azimuth: front face inner->outer, then back face outer->inner
    sec = []
    for i in range(n):
        pts = [sph(r, ef[k, i], az[i]) for k, r in enumerate(rhos)]
        pts += [sph(r, eb[k, i], az[i]) for k, r in reversed(list(enumerate(rhos)))]
        sec.append(pts)
    V = np.array(sec).reshape(-1, 3)
    L = 2 * m
    F = []
    for i in range(n):
        j = (i + 1) % n
        for e in range(L):
            a, b = i * L + e, i * L + (e + 1) % L
            c, d = j * L + (e + 1) % L, j * L + e
            F += [[a, b, c], [a, c, d]]
    mesh = trimesh.Trimesh(V, np.array(F), process=True)
    mesh.fix_normals()
    return mesh, (float(np.degrees(ef[-1].min())), float(np.degrees(eb[-1].max())))


def pin_el(band):
    """the pins sit at the middle of the ring's seat band"""
    return np.radians(sum(band) / 2)


def pin_meshes(band, r_extra=0.0):
    out = []
    el = pin_el(band)
    for az in (0.0, np.pi):
        u = sph(1.0, el, az)
        c = trimesh.creation.cylinder(radius=PIN_D / 2 + r_extra, height=PIN_IN + PIN_OUT, sections=32)
        c.apply_translation((0, 0, RHO_SEAT + (PIN_OUT - PIN_IN) / 2))
        c.apply_transform(to_matrix(rot_to(u)))
        out.append(c)
    return out


# ---------------------------------------------------------------- B-rep parts

def revolve_profile(pts):
    face = make_face(Polyline(*[(r, 0, z) for r, z in pts], close=True))
    return revolve(face, Axis.Z)


def sphere_arc(z0, z1, R, n=30):
    zs = np.linspace(z0, z1, n)
    return [(float(np.sqrt(max(R * R - z * z, 0))), float(z)) for z in zs]


def roller_sweep(row):
    """(max radius from the axis, z range) swept by a row's bearings as the carrier turns"""
    pts = []
    for uj in D.axes(row)[:1]:
        Rm = rot_to(uj)
        for a in np.linspace(0, 2 * np.pi, 48):
            for zz in (W.R - row.bearing.B / 2, W.R + row.bearing.B / 2):
                pts.append(Rm @ np.array([row.bearing.D / 2 * np.cos(a), row.bearing.D / 2 * np.sin(a), zz]))
    pts = np.array(pts)
    return float(np.hypot(pts[:, 0], pts[:, 1]).max()), float(pts[:, 2].min()), float(pts[:, 2].max())


def housings(band):
    """two halves split at z = 0: spherical seat, cavities for the rotating carriers, a floor with
    the input shaft bearing and a roof with the output carrier's main bearings"""
    r_in, z_in_lo, _ = roller_sweep(W.inp)
    r_ou, _, z_ou_hi = roller_sweep(W.out)
    R_seat = RHO_SEAT + SEAT_GAP
    cav_in, cav_out = r_in + 2.0, r_ou + 2.0
    z_seat_lo = -float(np.sqrt(max(R_seat ** 2 - cav_in ** 2, 0)))
    z_seat_hi = float(np.sqrt(max(R_seat ** 2 - cav_out ** 2, 0)))
    floor_in, roof = z_in_lo - 4.0, z_ou_hi + 4.0
    pts = [(R_OUT, 0.0), (R_OUT, floor_in - WALL), (SB["D"] / 2, floor_in - WALL), (SB["D"] / 2, floor_in),
           (cav_in, floor_in), (cav_in, z_seat_lo)] + sphere_arc(z_seat_lo, 0.0, R_seat)[1:]
    h_in = revolve_profile(pts)
    top = roof + WALL + OB["B"]
    pts = [(R_OUT, 0.0), (R_seat, 0.0)] + sphere_arc(0.0, z_seat_hi, R_seat)[1:] + \
          [(cav_out, roof), (OB["D"] / 2, roof), (OB["D"] / 2, top), (R_OUT, top)]
    h_out = revolve_profile(pts)
    return h_in, h_out, dict(floor_in=floor_in, roof=roof, top=top, z_seat_lo=z_seat_lo, z_seat_hi=z_seat_hi,
                             cav_in=cav_in, cav_out=cav_out)


def carrier(row, dims):
    """hub + a web along the row's cone + one stud per bearing, out through the bearing.
    Input: hub on the shaft. Output: hub tube through the main bearings, flange on top."""
    el = np.radians(abs(row.elev))
    c, s = np.cos(el), np.sin(el)
    h = WEB_T / 2
    # the web stops short of the ring's inner sphere (the back face comes down close to the output
    # cone); the studs bridge out to the bearings
    rw = W.R - row.bearing.B / 2 - 1.0 if row.elev < 0 else RHO_IN - 3.0
    if row.elev < 0:
        # web line (rho c, -rho s), normal (s, c)
        r0, r1 = SHAFT_D / 2, SHAFT_D / 2 + 8.0
        z_lo = dims["floor_in"] + 2.0
        q3 = (r1 + h * s) / c
        q6 = (r0 - h * s) / c
        pts = [(r0, z_lo), (r1, z_lo), (r1, -q3 * s - h * c), (rw * c - h * s, -rw * s - h * c),
               (rw * c + h * s, -rw * s + h * c), (r0, -q6 * s + h * c)]
    else:
        # web line (rho c, rho s), normal (-s, c)
        r0, r1 = OB["d"] / 2 - 4.0, OB["d"] / 2
        z_top = dims["top"] + 2.0
        q1 = (r0 - h * s) / c
        q4 = (r1 + h * s) / c
        pts = [(r0, q1 * s - h * c), (rw * c + h * s, rw * s - h * c), (rw * c - h * s, rw * s + h * c),
               (r1, q4 * s + h * c), (r1, z_top), (44.0, z_top), (44.0, z_top + 8.0), (r0, z_top + 8.0)]
    body = revolve_profile(pts)
    studs = []
    for uj in D.axes(row):
        L0, L1 = 26.0, W.R + row.bearing.B / 2 + 1.0
        studs.append(orient(Pos(0, 0, (L0 + L1) / 2) * Cylinder(row.bearing.d / 2, L1 - L0), uj))
    return body.fuse(*studs).clean()


def orient(shape, u):
    u = np.asarray(u, float)
    z = np.array([0, 0, 1.0])
    ax = np.cross(z, u)
    s = np.linalg.norm(ax)
    if s < 1e-12:
        return shape
    return shape.rotate(Axis((0, 0, 0), tuple(ax / s)), float(np.degrees(np.arctan2(s, z @ u))))


def bearing_shape(d, Dd, B):
    return Pos(0, 0, 0) * (Cylinder(Dd / 2, B) - Cylinder(d / 2, B + 1))


def roller_shape(row):
    """the bearing's outer ring as the roller: a tube D x B centred at R along +z"""
    return Pos(0, 0, W.R) * (Cylinder(row.bearing.D / 2, row.bearing.B) - Cylinder(row.bearing.d / 2 + 0.6, row.bearing.B + 1))


def input_shaft(dims):
    z0 = dims["floor_in"] - WALL - 16.0
    z1 = -18.0
    return Pos(0, 0, (z0 + z1) / 2) * Cylinder(SHAFT_D / 2, z1 - z0)


def motor(dims):
    z = dims["floor_in"] - WALL - 2.0 - MOTOR["L"] / 2
    return Pos(0, 0, z) * Cylinder(MOTOR["D"] / 2, MOTOR["L"])


# ---------------------------------------------------------------- meshes, poses, checks

def mesh_of(shape, tol=0.02, ang=0.1):
    OUT.mkdir(parents=True, exist_ok=True)
    tmp = OUT / "_t.stl"
    export_stl(shape, str(tmp), tolerance=tol, angular_tolerance=ang)
    m = trimesh.load(tmp)
    tmp.unlink()
    return m


def boolean(a, b, op):
    import manifold3d as mf
    def to_mf(m):
        return mf.Manifold(mf.Mesh(vert_properties=np.asarray(m.vertices, np.float32),
                                   tri_verts=np.asarray(m.faces, np.uint32)))
    A, B = to_mf(a), to_mf(b)
    r = {"diff": A - B, "and": A ^ B, "or": A + B}[op]
    o = r.to_mesh()
    return trimesh.Trimesh(np.asarray(o.vert_properties)[:, :3], np.asarray(o.tri_verts), process=False)


def union_all(ms):
    import manifold3d as mf
    out = None
    for m in ms:
        x = mf.Manifold(mf.Mesh(vert_properties=np.asarray(m.vertices, np.float32), tri_verts=np.asarray(m.faces, np.uint32)))
        out = x if out is None else out + x
    o = out.to_mesh()
    return trimesh.Trimesh(np.asarray(o.vert_properties)[:, :3], np.asarray(o.tri_verts), process=False)


def transformed(m, M):
    m = m.copy()
    m.apply_transform(M)
    return m


def pose(name, phi, u=None):
    """4x4 pose at input angle phi; parts are stored in their own frames (ring frame, carrier
    frames at angle 0)"""
    if name in ("housing_in", "housing_out", "motor", "shaft_bearing", "main_bearings"):
        return np.eye(4)
    if name == "ring":
        return to_matrix(D.ring_rot(W, phi))
    if name in ("input_carrier", "input_shaft"):
        return to_matrix(D.Rz(phi))
    if name == "output_carrier":
        return to_matrix(D.Rz(D.out_angle(W, phi)))
    if name == "input_roller":
        return to_matrix(D.Rz(phi) @ rot_to(u))
    if name == "output_roller":
        return to_matrix(D.Rz(D.out_angle(W, phi)) @ rot_to(u))
    raise ValueError(name)


def build():
    OUT.mkdir(parents=True, exist_ok=True)
    ring, band = ring_mesh()
    ring = union_all([ring] + pin_meshes(band))
    h_in, h_out, dims = housings(band)
    hin, hout = mesh_of(h_in), mesh_of(h_out)
    # pin slots: the pins' swept volume over a nutation cycle, with 0.1 mm clearance
    sweep = []
    for ps in np.linspace(0, 2 * np.pi, 72, endpoint=False):
        Rm = D.Rz(ps) @ D.Rx(np.radians(W.theta)) @ D.Rz(-ps)
        for p in pin_meshes(band, r_extra=0.1):
            sweep.append(transformed(p, to_matrix(Rm)))
    slot = union_all(sweep)
    hin, hout = boolean(hin, slot, "diff"), boolean(hout, slot, "diff")
    parts = dict(ring=ring, housing_in=hin, housing_out=hout,
                 input_carrier=mesh_of(carrier(W.inp, dims)), output_carrier=mesh_of(carrier(W.out, dims)),
                 input_roller=mesh_of(roller_shape(W.inp), tol=0.005, ang=0.05),
                 output_roller=mesh_of(roller_shape(W.out), tol=0.005, ang=0.05),
                 input_shaft=mesh_of(input_shaft(dims)), motor=mesh_of(motor(dims)),
                 shaft_bearing=mesh_of(Pos(0, 0, dims["floor_in"] - WALL / 2) * bearing_shape(SB["d"], SB["D"], SB["B"])),
                 main_bearings=mesh_of(Pos(0, 0, dims["roof"] + (WALL + OB["B"]) / 2) * bearing_shape(OB["d"], OB["D"], WALL + OB["B"])))
    return parts, dims, band


def rollers():
    """(name, carrier, u) for every bearing"""
    out = [("input_roller", uj) for uj in D.axes(W.inp)]
    out += [("output_roller", uj) for uj in D.axes(W.out)]
    return out


def check(parts, phis=(0.0, 0.4, 1.1, 2.3, 4.0)):
    res = dict(parts={})
    for n, m in parts.items():
        res["parts"][n] = dict(watertight=bool(m.is_watertight), bodies=len(m.split(only_watertight=False)),
                               volume_mm3=float(m.volume))
    fat = {}
    for n in ("input_roller", "output_roller"):
        f = parts[n].copy()
        f.vertices = f.vertices + f.vertex_normals * 0.03
        fat[n] = f
    worst_cut, untouched, n_roll = 0.0, 0, 0
    clash = {}
    for phi in phis:
        ring = transformed(parts["ring"], pose("ring", phi))
        for n, u in rollers():
            M = pose(n, phi, u)
            worst_cut = max(worst_cut, boolean(transformed(parts[n], M), ring, "and").volume)
            untouched += boolean(transformed(fat[n], M), ring, "and").volume < 1e-4
            n_roll += 1
        for other in ("housing_in", "housing_out", "input_carrier", "output_carrier", "input_shaft"):
            v = boolean(ring, transformed(parts[other], pose(other, phi)), "and").volume
            clash[f"ring/{other}"] = max(clash.get(f"ring/{other}", 0.0), v)
        for n, u in rollers():
            for hh in ("housing_in", "housing_out"):
                v = boolean(transformed(parts[n], pose(n, phi, u)), parts[hh], "and").volume
                clash[f"rollers/{hh}"] = max(clash.get(f"rollers/{hh}", 0.0), v)
        for c in ("input_carrier", "output_carrier"):
            for hh in ("housing_in", "housing_out"):
                v = boolean(transformed(parts[c], pose(c, phi)), parts[hh], "and").volume
                clash[f"{c}/{hh}"] = max(clash.get(f"{c}/{hh}", 0.0), v)
    res["interference"] = dict(roller_cut_mm3=worst_cut, rollers_not_touching=f"{untouched}/{n_roll}",
                               clashes_mm3={k: v for k, v in clash.items()})
    return res


def main():
    parts, dims, band = build()
    for n, m in parts.items():
        m.export(OUT / f"{n}.stl")
    res = check(parts)
    steel, alu = 7.85e-3, 2.70e-3
    mass = dict(ring=parts["ring"].volume * steel, housing=(parts["housing_in"].volume + parts["housing_out"].volume) * alu,
                input_carrier=parts["input_carrier"].volume * steel, output_carrier=parts["output_carrier"].volume * alu)
    mass["total (no motor, bought parts approx)"] = sum(mass.values()) + 3 * 19 + 12 * 12 + 35 + 2 * 25
    res["mass_g"] = mass
    res["dims"] = dims
    res["ring_band_deg"] = band
    res["bought"] = dict(input_rollers=f"3 x {W.inp.bearing.name}", output_rollers=f"12 x {W.out.bearing.name}",
                         shaft_bearing=SB, main_bearings=OB, pins=f"2 x dowel {PIN_D} mm", motor=MOTOR)
    (OUT / "parts.json").write_text(json.dumps(res, indent=1, default=float))
    for n, r in res["parts"].items():
        print(f"  {n:15s} watertight {r['watertight']!s:5s} bodies {r['bodies']}  {r['volume_mm3'] / 1e3:7.1f} cm3")
    print("  rollers:", res["interference"]["roller_cut_mm3"], "mm3 worst cut,", res["interference"]["rollers_not_touching"], "not touching")
    for k, v in res["interference"]["clashes_mm3"].items():
        print(f"  clash {k:32s} {v:.4f} mm3")
    print("  mass (g):", {k: round(v) for k, v in mass.items()})
    bad = [n for n, r in res["parts"].items() if not r["watertight"] or r["bodies"] != 1]
    print("PART CHECK", "FAILED " + str(bad) if bad else "ok")


if __name__ == "__main__":
    main()
