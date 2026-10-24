"""
Presentation renders of the roller-screw actuator: cad.py's assembly -> GLB with PBR
materials -> Blender (Cycles), using the product-design skill's blender_render.py.

Scenes
  section   half-section through the axis, lying along x: rotor, nut and rollers inside
            the motor; the screw and the rollers behind the plane left whole
  screw     the roller screw alone, the nut and its end parts cut in half

Viz only: materials and cuts live here, not in cad.py.

    python rollerscrew/render.py        (or: make rollerscrew-render)
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import trimesh
from build123d import Align, Box, Compound, Pos, Rot, export_stl
from trimesh.visual.material import PBRMaterial

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cad  # noqa: E402

OUT = HERE / "out" / "render"
BLENDER = os.environ.get("BLENDER", "/opt/blender-5.0.1-linux-x64/blender")
SCRIPT = Path.home() / ".claude/skills/product-design/scripts/blender_render.py"

STEEL = ((0.62, 0.64, 0.67, 1), 1.0, 0.28)
DARK = ((0.10, 0.11, 0.13, 1), 0.6, 0.45)
MATERIALS = {
    "screw": ((0.72, 0.74, 0.77, 1), 1.0, 0.22), "rod_eye": STEEL,
    "roller": ((0.80, 0.70, 0.45, 1), 1.0, 0.25), "ring_gear": ((0.70, 0.62, 0.40, 1), 1.0, 0.3),
    "carrier": ((0.50, 0.52, 0.55, 1), 1.0, 0.35), "nut": ((0.58, 0.52, 0.42, 1), 1.0, 0.3),
    "rotor": ((0.16, 0.17, 0.19, 1), 0.9, 0.4), "magnet": ((0.78, 0.80, 0.82, 1), 1.0, 0.18),
    "stator_core": ((0.22, 0.24, 0.27, 1), 0.9, 0.45), "coil": ((0.85, 0.45, 0.22, 1), 1.0, 0.28),
    "bearing_front": ((0.85, 0.86, 0.88, 1), 1.0, 0.15), "bearing_rear": ((0.85, 0.86, 0.88, 1), 1.0, 0.15),
    "housing": ((0.70, 0.72, 0.75, 1), 1.0, 0.38), "cap_front": DARK, "cap_rear": DARK, "tail_tube": DARK,
    "load_cell": ((0.80, 0.62, 0.18, 1), 1.0, 0.3), "clevis": STEEL,
    "encoder_magnet": ((0.30, 0.45, 0.70, 1), 0.3, 0.4), "encoder_board": ((0.12, 0.40, 0.28, 1), 0.3, 0.5),
    "boot": ((0.05, 0.05, 0.05, 1), 0.0, 0.8),
}


def solids_of(r):
    parts = list(r) if isinstance(r, (list, tuple)) else [r]
    return [x for p in parts for x in (p.solids() if hasattr(p, "solids") else [p])]


def cut(solids):
    """keep the y >= 0 half (the camera looks from -y)"""
    keep = Pos(0, 0, -50) * Box(200, 100, 400, align=(Align.CENTER, Align.MIN, Align.MIN))
    out = []
    for s in solids:
        bb = s.bounding_box()
        if bb.min.Y >= 0:
            out.append(s)
        elif bb.max.Y > 0:
            out += solids_of(s & keep)
    return out


LAY = Rot(0, 90, 0)          # actuator axis z -> x


def scene_section():
    p = cad.assembly()
    s = {}
    whole = ("screw", "rod_eye", "clevis")
    for k, v in p.items():
        if k in whole:
            s[k] = v
        elif k == "roller":
            s[k] = [r for r in v if r.bounding_box().center().Y > -0.5]
        else:
            s[k] = cut(v)
    return {k: [LAY * x for x in v] for k, v in s.items() if v}


def scene_screw():
    p = cad.assembly()
    s = {k: p[k] for k in ("screw", "roller")}
    for k in ("nut", "ring_gear", "carrier"):
        s[k] = cut(p[k])
    z = (cad.Z_NUT0 + cad.Z_NUT1) / 2
    clip = Pos(0, 0, z - 40) * Box(60, 60, 80, align=(Align.CENTER, Align.CENTER, Align.MIN))
    s["screw"] = solids_of(s["screw"][0] & clip)
    return {k: [LAY * (Pos(0, 0, -z) * x) for x in v] for k, v in s.items()}


def to_glb(name, parts):
    OUT.mkdir(parents=True, exist_ok=True)
    scene = trimesh.Scene()
    for part, solids in parts.items():
        stl = OUT / f"_{name}_{part}.stl"
        export_stl(Compound(children=list(solids)), str(stl), tolerance=0.01, angular_tolerance=0.1)
        mesh = trimesh.load(stl)
        rgba, metal, rough = MATERIALS[part]
        mesh.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=part, baseColorFactor=rgba, metallicFactor=metal, roughnessFactor=rough,
            alphaMode="OPAQUE", doubleSided=True))
        scene.add_geometry(mesh, geom_name=part)
        stl.unlink()
    scene.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))   # Z-up -> Y-up
    glb = OUT / f"{name}.glb"
    scene.export(glb)
    return glb


def blender(glb, png, samples=160, res=1500):
    if png.exists():
        png.unlink()
    env = dict(os.environ, LC_ALL="C", LANG="C")
    subprocess.run([BLENDER, "-b", "-P", str(SCRIPT), "--", str(glb), str(png), "--samples", str(samples),
                    "--res", str(res), "--bg"], env=env, check=True, capture_output=True)
    if not png.exists():
        raise RuntimeError(f"Blender wrote nothing for {png}")
    print(f"wrote {png}")


def main():
    for name, build in (("section", scene_section), ("screw", scene_screw)):
        glb = to_glb(name, build())
        print(f"{glb.name}: {glb.stat().st_size / 1e6:.1f} MB")
        blender(glb, OUT / f"{name}.png")


if __name__ == "__main__":
    main()
