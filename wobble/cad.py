"""
CAD of the twin wobble reducer (design.py), placed by its own kinematics.

Layout along the axis (z), a clevis joint:

    motor | end housing A (fixed ring 1) | output ring (both output grooves) | end housing B (fixed ring 2)
                     block 1 (z = -20)                    block 2 (z = +20, mirrored)

The two end housings are tied by a U frame (the clevis); the output ring sits between them with
its lug out of the open side of the U and NO bearing of its own: the cam rollers locate it.

Every surface a roller touches comes from design.py: the groove walls are the spherical offsets
of the roller paths, extruded radially through the block centre (rollers are cones with apex
there). Those walls are cut with a mesh boolean (manifold3d), the build123d-part fallback for
geometry B-rep won't merge cleanly: a periodic ruled wave needs ~3000 stations to hold the
1 mm crest radius. Everything else is B-rep.

    python wobble/cad.py            (or: make wobble-cad)   -> out/cad/*.stl, parts.json, checks
"""

import json
import sys
from pathlib import Path

import numpy as np
import trimesh
from build123d import (Axis, Box, Cylinder, Location, Plane, Polyline, Pos, Rot, Solid, export_stl,
                       make_face, revolve, Compound)

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "cad"
sys.path.insert(0, str(HERE))
import design as D  # noqa: E402

W = D.Wobble()

# ---------------------------------------------------------------- dimensions (mm)
RHO0, RHO1 = 34.0, 46.0          # roller cone, distance from the block centre (pitch R = 40 in the middle)
STEM_D, STEM_IN = 3.0, 24.0      # roller stem (runs in a drilled hole in the block) and its inner end
BLOCK_SR = 33.0                  # block outer sphere
BLOCK_H = 13.5                   # block half height (truncated sphere)
GAP_SPH = 0.6                    # block sphere to ring bore
RING_SR = BLOCK_SR + GAP_SPH     # ring bore sphere, concentric with the block: constant gap while it nutates
GROOVE_DEEP = RHO1 + 0.8         # groove bottom (roller end clearance)
RING_R = 58.0                    # housing / output ring outer radius
ZC = abs(W.blocks[0][0])         # block centres at -ZC, +ZC
WALL_IN, WALL_T = 36.5, 6.0      # end wall inner face (|z|) and thickness
WB = dict(name="6806-2RS thin-section ball bearing", d=30.0, D=42.0, B=7.0)   # wobble bearing, one per block
SB = dict(name="6802-2RS ball bearing", d=15.0, D=24.0, B=5.0)                # input shaft, one per end wall
LUG = (70.0, 30.0, 16.0)         # output lug: reach beyond the ring, width, thickness
FRAME_T, FRAME_GAP = 8.0, 4.0    # clevis plates and the gap to the output ring
MOTOR = dict(name="pancake BLDC (stand-in)", D=72.0, L=26.0)


def block_place(k):
    """(A, o) of block k: the canonical block mapped by A then moved to o (design.placements)"""
    return D.placements(W)[k]


def to_matrix(A, o, R=np.eye(3)):
    M = np.eye(4)
    M[:3, :3] = A @ R
    M[:3, 3] = o
    return M


# ---------------------------------------------------------------- grooves (mesh)

def wall_curves(row, n=3600):
    """the groove's lower (A) and upper (B) walls on the unit sphere of the ring frame, as
    elevation (rad) on a uniform azimuth grid"""
    gv = D.groove(W, row, samples=12000)
    az = np.linspace(0, 2 * np.pi, n, endpoint=False)
    out = []
    for key in ("A", "B"):
        a = np.mod(np.radians(gv[key]["az"]), 2 * np.pi)
        e = np.radians(gv[key]["el"])
        o = np.argsort(a)
        a, e = a[o], e[o]
        a = np.concatenate([a[-1:] - 2 * np.pi, a, a[:1] + 2 * np.pi])
        e = np.concatenate([e[-1:], e, e[:1]])
        out.append(np.interp(az, a, e))
    return az, out[0], out[1]


def sph(rho, el, az):
    return np.stack([rho * np.cos(el) * np.cos(az), rho * np.cos(el) * np.sin(az), rho * np.sin(el)], -1)


def groove_mesh(row, rho0=RING_SR - 1.0, rho1=GROOVE_DEEP):
    """closed tube: between the two walls in elevation, rho0..rho1 radially (block-centred frame)"""
    az, lo, hi = wall_curves(row)
    n = len(az)
    v = np.concatenate([sph(rho0, lo, az), sph(rho1, lo, az), sph(rho1, hi, az), sph(rho0, hi, az)])
    f = []
    for i in range(n):
        j = (i + 1) % n
        for s in range(4):
            a, b = s * n, ((s + 1) % 4) * n
            f += [[a + i, b + i, b + j], [a + i, b + j, a + j]]
    m = trimesh.Trimesh(v, np.array(f), process=True)
    m.fix_normals()
    return m


# ---------------------------------------------------------------- B-rep parts

def revolve_profile(pts):
    """revolve an (r, z) polygon about z"""
    face = make_face(Polyline(*[(r, 0, z) for r, z in pts], close=True))
    return revolve(face, Axis.Z)


def sphere_arc(zc, z0, z1, R, n=24):
    zs = np.linspace(z0, z1, n)
    return [(float(np.sqrt(max(R * R - (z - zc) ** 2, 0))), float(z)) for z in zs]


def end_housing(side):
    """fixed ring + end wall + input bearing seat; side -1 (motor end, block 1) or +1 (block 2).
    Built for side -1 and mirrored."""
    zc = -ZC
    z_end = -(WALL_IN + WALL_T)
    # the bore sphere meets the end wall face at radius sqrt(R^2 - (wall - zc)^2)
    rw = float(np.sqrt(RING_SR ** 2 - (WALL_IN - ZC) ** 2))
    pts = [(RING_R, zc), (RING_R, z_end), (SB["D"] / 2, z_end), (SB["D"] / 2, -WALL_IN), (rw, -WALL_IN)]
    pts += sphere_arc(zc, -WALL_IN, zc, RING_SR)[1:]
    p = revolve_profile(pts)
    # bolt holes to the clevis plate
    for a in range(6):
        ang = np.radians(30 + 60 * a)
        p -= Pos(50 * np.cos(ang), 50 * np.sin(ang), z_end + 3) * Cylinder(1.7, 10)
    return p if side < 0 else p.mirror(Plane.XY)


def output_ring():
    """middle ring: two spherical bores (one per block), the lug"""
    pts = [(RING_R, -ZC), (RING_R, ZC)]
    pts += sphere_arc(ZC, ZC, 0.0, RING_SR)[1:]
    pts += sphere_arc(-ZC, 0.0, -ZC, RING_SR)[1:]
    ring = revolve_profile(pts)
    reach, wid, thk = LUG
    lug = Pos(0, RING_R - 3 + reach / 2, 0) * Box(wid, reach, thk)
    lug -= Pos(0, RING_R + reach - 12, 0) * Cylinder(4.2, thk + 2)
    return ring.fuse(lug).clean()


def block():
    """canonical block (centre at origin): truncated sphere, wobble-bearing bore, stem holes"""
    b = Solid.make_sphere(BLOCK_SR) & Box(2 * BLOCK_SR + 2, 2 * BLOCK_SR + 2, 2 * BLOCK_H)
    b -= Cylinder(WB["D"] / 2, 2 * BLOCK_H + 2)
    holes = []
    for row in (W.fixed, W.output):
        u, _ = D.axes(W, row)
        for uj in u:
            holes.append(axis_cylinder(uj, STEM_D / 2 + 0.05, WB["D"] / 2 - 1, BLOCK_SR + 1))
    return b.cut(*holes)


def axis_cylinder(u, r, d0, d1):
    """cylinder of radius r along unit vector u from distance d0 to d1 from the origin"""
    c = Pos(0, 0, (d0 + d1) / 2) * Cylinder(r, d1 - d0)
    return _orient(c, u)


def _orient(shape, u):
    """rotate a shape built along +z so +z points along u"""
    u = np.asarray(u, float)
    z = np.array([0, 0, 1.0])
    ax = np.cross(z, u)
    s = np.linalg.norm(ax)
    if s < 1e-12:
        return shape
    ang = np.degrees(np.arctan2(s, z @ u))
    return shape.rotate(Axis((0, 0, 0), tuple(ax / s)), ang)


def roller():
    """cone frustum (apex at the block centre) plus its stem, along +z"""
    g = W.gamma
    d0, d1 = RHO0 * np.cos(g), RHO1 * np.cos(g)
    cone = Pos(0, 0, d0) * Solid.make_cone(d0 * np.tan(g), d1 * np.tan(g), d1 - d0)
    stem = Pos(0, 0, (STEM_IN + d0 + 0.5) / 2) * Cylinder(STEM_D / 2, d0 + 0.5 - STEM_IN)
    return cone.fuse(stem).clean()


def bearing(spec):
    return Cylinder(spec["D"] / 2, spec["B"]) - Cylinder(spec["d"] / 2, spec["B"] + 1)


def input_shaft():
    """hollow-free shaft, z from motor coupling to the far end wall, with two tilted journals.
    Journal k carries the wobble bearing of block k, tilted with that block at phi = 0."""
    z0, z1 = -(WALL_IN + WALL_T + 14), WALL_IN + WALL_T - 0.5
    core = Pos(0, 0, (z0 + z1) / 2) * Cylinder(SB["d"] / 2, z1 - z0)
    mid = Pos(0, 0, 0) * Cylinder(10.0, 2 * (ZC + 6))
    parts = [core, mid]
    for k in range(2):
        A, o = block_place(k)
        Rb = D.block_rot(W, 0.0)
        # the journal is the block's bearing seat: a cylinder about the block axis, built in the block frame
        j = Cylinder(WB["d"] / 2, WB["B"] + 2)
        parts.append(place_shape(j, A @ Rb, o))
    return parts[0].fuse(*parts[1:]).clean()


def place_shape(shape, Rm, o):
    """apply rotation/reflection matrix Rm then translation o to a build123d shape"""
    from OCP.gp import gp_GTrsf, gp_Mat, gp_XYZ
    from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform, BRepBuilderAPI_Transform
    from OCP.gp import gp_Trsf
    det = np.linalg.det(Rm)
    if det > 0:
        t = gp_Trsf()
        t.SetValues(*Rm[0], o[0], *Rm[1], o[1], *Rm[2], o[2])
        return Solid(BRepBuilderAPI_Transform(shape.wrapped, t, True).Shape())
    # a reflection: mirror in z first, then the remaining rotation
    Mz = np.diag([1, 1, -1.0])
    s = shape.mirror(Plane.XY)
    return place_shape(s, Rm @ Mz, o)


def clevis():
    """U frame: a plate on each end housing, joined by a bridge on the -y side"""
    zo = WALL_IN + WALL_T
    plate_r = RING_R + 2
    bridge_y = -(RING_R + FRAME_GAP + 8)
    parts = []
    for s in (-1, 1):
        z = s * (zo + FRAME_T / 2)
        disc = Pos(0, 0, z) * Cylinder(plate_r, FRAME_T)
        arm = Pos(0, (bridge_y - 0) / 2, z) * Box(2 * plate_r * 0.8, abs(bridge_y), FRAME_T)
        p = disc.fuse(arm)
        p -= Pos(0, 0, z) * Cylinder(SB["D"] / 2 + 2, FRAME_T + 2)
        parts.append(p)
    bridge = Pos(0, bridge_y, 0) * Box(2 * plate_r * 0.8, 16, 2 * (zo + FRAME_T))
    return parts[0].fuse(parts[1], bridge).clean()


def motor():
    zo = WALL_IN + WALL_T + FRAME_T
    return Pos(0, 0, -(zo + 2 + MOTOR["L"] / 2)) * Cylinder(MOTOR["D"] / 2, MOTOR["L"])


# ---------------------------------------------------------------- meshes

def mesh_of(shape, tol=0.02, ang=0.1):
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
    out = r.to_mesh()
    return trimesh.Trimesh(np.asarray(out.vert_properties)[:, :3], np.asarray(out.tri_verts), process=False)


def transformed(m, M):
    m = m.copy()
    m.apply_transform(M)              # trimesh flips the winding itself for a reflection
    return m


def build():
    """every part as a mesh in its assembled pose at phi = 0, plus the pieces the sim needs"""
    OUT.mkdir(parents=True, exist_ok=True)
    gm = {"fixed": groove_mesh(W.fixed), "output": groove_mesh(W.output)}
    parts = {}
    # fixed rings in the end housings
    for k, side in ((0, -1), (1, 1)):
        A, o = block_place(k)
        h = mesh_of(end_housing(side))
        h = boolean(h, transformed(gm["fixed"], to_matrix(A, o)), "diff")
        parts[f"housing_{'ab'[k]}"] = h
    ring = mesh_of(output_ring())
    for k in (0, 1):
        A, o = block_place(k)
        ring = boolean(ring, transformed(gm["output"], to_matrix(A, o)), "diff")
    parts["output_ring"] = ring
    blk = mesh_of(block())
    rol = mesh_of(roller(), tol=0.005, ang=0.05)
    parts["block_canonical"] = blk
    parts["roller_canonical"] = rol
    parts["input_shaft"] = mesh_of(input_shaft())
    parts["wobble_bearing"] = mesh_of(bearing(WB))
    parts["shaft_bearing"] = mesh_of(bearing(SB))
    parts["clevis"] = mesh_of(clevis())
    parts["motor"] = mesh_of(motor())
    return parts


def roller_axes():
    """(row, block index, canonical axis u) for every roller"""
    out = []
    for k in range(len(W.blocks)):
        for row, nm in ((W.fixed, "fixed"), (W.output, "output")):
            u, side = D.axes(W, row)
            for j, uj in enumerate(u):
                out.append((nm, k, j, uj, int(side[j])))
    return out


def rot_to(u):
    """rotation matrix taking +z to u"""
    z = np.array([0, 0, 1.0])
    ax = np.cross(z, u)
    s = np.linalg.norm(ax)
    if s < 1e-12:
        return np.eye(3)
    return trimesh.transformations.rotation_matrix(np.arctan2(s, z @ u), ax / s)[:3, :3]


def pose(name, phi, k=None, u=None):
    """4x4 pose of a part at input angle phi (assembled frame = phi 0 for rings/housings)"""
    if name in ("housing_a", "housing_b", "clevis", "motor"):
        return np.eye(4)
    if name == "output_ring":
        M = np.eye(4)
        M[:3, :3] = D.Rz(D.output_angle(W, phi))
        return M
    if name == "input_shaft":
        M = np.eye(4)
        M[:3, :3] = D.Rz(phi)
        return M
    A, o = block_place(k)
    Rb = D.block_rot(W, phi)
    if name == "block":
        return to_matrix(A, o, Rb)
    if name == "roller":
        return to_matrix(A, o, Rb @ rot_to(u))
    if name == "wobble_bearing":
        return to_matrix(A, o, Rb)
    raise ValueError(name)


def check(parts, phis=(0.0, 0.7, 1.9, 3.3, 5.1)):
    """part checks, then interference: rollers vs grooves (must touch, must not cut), blocks vs
    rings (must clear), at several input angles"""
    import manifold3d as mf  # noqa: F401
    res = dict(parts={}, interference=[])
    for n, m in parts.items():
        res["parts"][n] = dict(watertight=bool(m.is_watertight), bodies=len(m.split(only_watertight=False)),
                               volume_mm3=float(m.volume))
    rings = {"housing_a": parts["housing_a"], "housing_b": parts["housing_b"], "output_ring": parts["output_ring"]}
    ring_of = {("fixed", 0): "housing_a", ("fixed", 1): "housing_b", ("output", 0): "output_ring", ("output", 1): "output_ring"}
    rol = parts["roller_canonical"]
    rol_fat = rol.copy()
    rol_fat.vertices = rol_fat.vertices + rol_fat.vertex_normals * 0.03          # 30 um proud
    worst_cut, untouched, blk_hit = 0.0, 0, 0.0
    for phi in phis:
        ringpose = {n: transformed(m, pose(n, phi)) for n, m in rings.items()}
        for row, k, j, u, side in roller_axes():
            Mr = pose("roller", phi, k, u)
            ring = ringpose[ring_of[(row, k)]]
            cut = boolean(transformed(rol, Mr), ring, "and").volume
            touch = boolean(transformed(rol_fat, Mr), ring, "and").volume
            worst_cut = max(worst_cut, cut)
            untouched += touch < 1e-4
        for k in (0, 1):
            Mb = pose("block", phi, k)
            b = transformed(parts["block_canonical"], Mb)
            for n, ring in ringpose.items():
                blk_hit = max(blk_hit, boolean(b, ring, "and").volume)
    n_roll = len(roller_axes()) * len(phis)
    res["interference"] = dict(roller_cut_mm3=worst_cut, rollers_not_touching=f"{untouched}/{n_roll}",
                               block_vs_rings_mm3=blk_hit)
    return res


def main():
    parts = build()
    for n, m in parts.items():
        m.export(OUT / f"{n}.stl")
    res = check(parts)
    steel, alu = 7.85e-3, 2.70e-3                                         # g/mm3
    mats = dict(housing_a=alu, housing_b=alu, output_ring=steel, block_canonical=steel, roller_canonical=steel,
                input_shaft=steel, clevis=alu)
    mass = {n: parts[n].volume * mats[n] for n in mats}
    mass["total (2 blocks, 44 rollers)"] = (mass["housing_a"] + mass["housing_b"] + mass["output_ring"] + mass["clevis"]
                                            + 2 * mass["block_canonical"] + 44 * mass["roller_canonical"]
                                            + mass["input_shaft"])
    res["mass_g"] = mass
    res["bought"] = dict(wobble_bearing=WB, shaft_bearing=SB, motor=MOTOR)
    (OUT / "parts.json").write_text(json.dumps(res, indent=1, default=float))
    for n, r in res["parts"].items():
        print(f"  {n:18s} watertight {r['watertight']!s:5s} bodies {r['bodies']}  {r['volume_mm3'] / 1e3:7.1f} cm3")
    print("  interference:", res["interference"])
    print("  mass (g):", {k: round(v) for k, v in mass.items()})
    bad = [n for n, r in res["parts"].items() if not r["watertight"] or r["bodies"] != 1]
    print("PART CHECK", "FAILED " + str(bad) if bad else "ok")


if __name__ == "__main__":
    main()
