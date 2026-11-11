"""
Meshes and motion data for the interactive assembly page (web/index.html): cad.py's parts at a
web-friendly resolution, an exploded stacking along the axis, and the reducer's kinematics
(design.py: input carrier Rz(phi), ring Rz(psi) Rx(theta) Rz(-psi), output Rz(phi / i + clock))
posed live in the browser.

    python wobble2/web_export.py      (or: make wobble2-web)  -> web/parts.gltf.json, web/meta.json

then serve the folder (python -m http.server -d wobble2/web) and open http://localhost:8000
"""

import json
import sys
from pathlib import Path

import numpy as np
import trimesh
from trimesh.visual.material import PBRMaterial

HERE = Path(__file__).resolve().parent
WEB = HERE / "web"
sys.path.insert(0, str(HERE))
import cad as C  # noqa: E402
import design as D  # noqa: E402

W = C.W
GAP = 18.0                     # mm between exploded groups along the axis

# colour (sRGB), metallic, roughness: render.py's palette
ALU = ((0.80, 0.81, 0.83, 1), 0.9, 0.35)
STEEL = ((0.62, 0.64, 0.67, 1), 1.0, 0.25)
DARK = ((0.16, 0.16, 0.18, 1), 0.5, 0.45)
ORANGE = ((0.93, 0.45, 0.12, 1), 0.0, 0.45)
BLUE = ((0.20, 0.40, 0.70, 1), 0.3, 0.35)
BRG = ((0.82, 0.84, 0.87, 1), 1.0, 0.2)

# assembly order, bottom (motor) to top (output flange)
GROUPS = [
    # key, name, description, qty, kind, colour, fade (made see-through by "See inside")
    ("motor", "motor", f"{C.MOTOR['name']}, Ø{C.MOTOR['D']:.0f}×{C.MOTOR['L']:.0f}", 1, "bought", DARK, False),
    ("housing_in", "housing, input half", "lower half of the spherical seat, shaft bearing bore, pin slots", 1, "machined", ALU, True),
    ("shaft_bearing", "input shaft bearing", f"{C.SB['name']}, {C.SB['d']:.0f}×{C.SB['D']:.0f}×{C.SB['B']:.0f}", 1, "bought", BRG, False),
    ("input_shaft", "input shaft", f"Ø{C.SHAFT_D:.0f} shaft from the motor", 1, "machined", STEEL, False),
    ("input_carrier", "input carrier", f"hub, web and {len(W.inp.idx)} studs on a {W.inp.elev:+.0f}° cone; turns with the input",
     1, "machined", STEEL, False),
    ("input_rollers", "input rollers", f"{W.inp.bearing.name} ball bearings rolling on the {W.inp.lobes}-lobe front wave",
     len(W.inp.idx), "bought", BRG, False),
    ("ring", "wobble ring", f"band of a sphere, {W.inp.lobes}-lobe front wave, {W.out.lobes}-lobe back wave, 2 anti-rotation "
     f"pins; nutates {W.theta}° and never spins", 1, "machined", ORANGE, False),
    ("output_rollers", "output rollers", f"{W.out.bearing.name} on the {W.out.lobes}-lobe back wave", len(W.out.idx), "bought", BRG, False),
    ("output_carrier", "output carrier", f"{len(W.out.idx)} studs on a {W.out.elev:+.0f}° cone, hub tube, output flange", 1, "machined", BLUE, False),
    ("housing_out", "housing, output half", "upper half of the seat, roof with the main-bearing bore, pin slots", 1, "machined", ALU, True),
    ("main_bearings", "output main bearings", f"{C.OB['name']}, {C.OB['d']:.0f}×{C.OB['D']:.0f}×{C.OB['B']:.0f}", 2, "bought", BRG, False),
]
ROLE = dict(housing_in="fixed", housing_out="fixed", motor="fixed", shaft_bearing="fixed", main_bearings="fixed",
            input_shaft="input", input_carrier="input", input_rollers="input", ring="ring",
            output_rollers="output", output_carrier="output")


def web_parts():
    """cad.build() at web resolution: 720 azimuth stations on the ring faces (2880 for the
    checks), coarser B-rep tessellation, the pin slots swept over 24 steps (72)"""
    C.face_curves.__defaults__ = (720, 9000)
    ring, band = C.ring_mesh()
    ring = C.union_all([ring] + C.pin_meshes(band))
    h_in, h_out, dims = C.housings(band)
    coarse = lambda s, t=0.05: C.mesh_of(s, tol=t, ang=0.3)  # noqa: E731
    sweep = []
    for ps in np.linspace(0, 2 * np.pi, 24, endpoint=False):
        Rm = D.Rz(ps) @ D.Rx(np.radians(W.theta)) @ D.Rz(-ps)
        sweep += [C.transformed(p, C.to_matrix(Rm)) for p in C.pin_meshes(band, r_extra=0.15)]
    slot = C.union_all(sweep)
    return dict(
        ring=ring,
        housing_in=C.boolean(coarse(h_in), slot, "diff"), housing_out=C.boolean(coarse(h_out), slot, "diff"),
        input_carrier=coarse(C.carrier(W.inp, dims)), output_carrier=coarse(C.carrier(W.out, dims)),
        input_rollers=coarse(C.roller_shape(W.inp), 0.02), output_rollers=coarse(C.roller_shape(W.out), 0.02),
        input_shaft=coarse(C.input_shaft(dims)), motor=coarse(C.motor(dims)),
        shaft_bearing=coarse(C.Pos(0, 0, dims["floor_in"] - C.WALL / 2) * C.bearing_shape(C.SB["d"], C.SB["D"], C.SB["B"])),
        main_bearings=coarse(C.Pos(0, 0, dims["roof"] + (C.WALL + C.OB["B"]) / 2)
                             * C.bearing_shape(C.OB["d"], C.OB["D"], C.WALL + C.OB["B"])))


def objects(meshes):
    """(name, group, role, base 4x4, mesh): one per part; a roller row is one object per bearing,
    its base the turn from +z onto that bearing's axis"""
    rows = []
    for g, (key, *_r) in enumerate(GROUPS):
        if key in ("input_rollers", "output_rollers"):
            row = W.inp if key == "input_rollers" else W.out
            for j, u in enumerate(D.axes(row)):
                rows.append((f"{key}_{j}", g, ROLE[key], C.to_matrix(C.rot_to(u)), meshes[key]))
        else:
            rows.append((key, g, ROLE[key], np.eye(4), meshes[key]))
    return rows


def pose(role, base, phi):
    R = dict(fixed=np.eye(3), input=D.Rz(phi), ring=D.ring_rot(W, phi), output=D.Rz(D.out_angle(W, phi)))[role]
    return C.to_matrix(R) @ base


def explode_offsets(rows):
    """each group starts GAP above where the previous one ends (posed at phi = 0)"""
    dz, z_next = {}, None
    for g in range(len(GROUPS)):
        b = np.array([C.transformed(m, pose(role, base, 0.0)).bounds for _n, gg, role, base, m in rows if gg == g])
        lo, hi = b[:, 0, 2].min(), b[:, 1, 2].max()
        dz[g] = 0.0 if z_next is None else max(0.0, z_next - lo)
        z_next = hi + dz[g] + GAP
    return dz


def to_gltf(rows, path):
    scene = trimesh.Scene()
    for name, g, _role, _base, m in rows:
        col, metal, rough = GROUPS[g][5]
        m = m.copy()
        m.merge_vertices()
        m.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=name, baseColorFactor=col, metallicFactor=metal, roughnessFactor=rough, doubleSided=True))
        scene.add_geometry(m, geom_name=name, node_name=name)            # part frame; the page poses it
    from trimesh.exchange.gltf import export_gltf
    path.write_bytes(export_gltf(scene, embed_buffers=True)["model.gltf"])


def main():
    WEB.mkdir(exist_ok=True)
    rows = objects(web_parts())
    dz = explode_offsets(rows)
    to_gltf(rows, WEB / "parts.gltf.json")
    meta = dict(
        ratio=W.ratio, theta=W.theta, psi_k=W.inp.s * W.inp.slots, back_clock=W.back_clock, R=W.R,
        inp=dict(n=len(W.inp.idx), lobes=W.inp.lobes, elev=W.inp.elev, bearing=W.inp.bearing.name),
        out=dict(n=len(W.out.idx), lobes=W.out.lobes, elev=W.out.elev, bearing=W.out.bearing.name),
        groups=[dict(key=k, name=n, text=t, qty=q, kind=kd, dz=dz[g]) for g, (k, n, t, q, kd, _c, _f) in enumerate(GROUPS)],
        objects=[dict(name=n, group=g, role=role, fade=GROUPS[g][6], base=base.T.flatten().round(9).tolist())
                 for n, g, role, base, _m in rows])
    (WEB / "meta.json").write_text(json.dumps(meta, indent=1))
    size = (WEB / "parts.gltf.json").stat().st_size / 1e6
    tris = sum(len(r[4].faces) for r in rows)
    print(f"wrote {WEB}/parts.gltf.json ({len(rows)} parts, {tris} triangles, {size:.1f} MB) and meta.json "
          f"({len(GROUPS)} assembly steps, ratio {W.ratio:g})")


if __name__ == "__main__":
    main()
