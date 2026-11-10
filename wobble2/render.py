"""
Renders of wobble2 from cad.py's meshes (product-design pipeline: world-placed meshes -> PBR GLB
-> headless Blender with the product-design skill's studio script).

    python wobble2/render.py [--quick]   (or: make wobble2-render)
        -> out/render_full.png, out/render_cut.png, out/render_core.png (the ring and its two roller sets)
"""

import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import trimesh
from trimesh.visual.material import PBRMaterial

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cad as C  # noqa: E402

BLENDER = os.environ.get("BLENDER", "/opt/blender-5.0.1-linux-x64/blender")
SCRIPT = Path.home() / ".claude/skills/product-design/scripts/blender_render.py"
CAD = HERE / "out" / "cad"

ALU = ((0.80, 0.81, 0.83, 1), 0.9, 0.35)
STEEL = ((0.62, 0.64, 0.67, 1), 1.0, 0.25)
DARK = ((0.12, 0.12, 0.13, 1), 0.5, 0.45)
ORANGE = ((0.93, 0.45, 0.12, 1), 0.0, 0.45)
BLUE = ((0.20, 0.40, 0.70, 1), 0.3, 0.35)
BRG = ((0.82, 0.84, 0.87, 1), 1.0, 0.2)
CUT = ("housing_in", "housing_out")


def scene_items(phi=0.0):
    m = {p.stem: trimesh.load(p) for p in CAD.glob("*.stl")}
    items = [("housing_in", ALU), ("housing_out", ALU), ("ring", ORANGE), ("input_carrier", STEEL),
             ("output_carrier", BLUE), ("input_shaft", STEEL), ("motor", DARK), ("shaft_bearing", BRG),
             ("main_bearings", BRG)]
    out = [(n, C.transformed(m[n], C.pose(n, phi)), mat) for n, mat in items]
    rol = [C.transformed(m[n], C.pose(n, phi, u)) for n, u in C.rollers()]
    out.append(("rollers", trimesh.util.concatenate(rol), BRG))
    return out


def cutaway(items):
    box = trimesh.creation.box(extents=(400, 400, 400))
    box.apply_translation((0, 200, 0))
    return [(n, C.boolean(mesh, box, "diff") if n in CUT else mesh, mat) for n, mesh, mat in items]


def glb(items, path):
    sc = trimesh.Scene()
    for n, mesh, (col, metal, rough) in items:
        mesh = mesh.copy()
        mesh.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=n, baseColorFactor=col, metallicFactor=metal, roughnessFactor=rough, doubleSided=True))
        sc.add_geometry(mesh, geom_name=n, node_name=n)
    # cut side toward the camera (it looks from -y), axis horizontal, turned a little, Z-up -> Y-up
    sc.apply_transform(trimesh.transformations.rotation_matrix(np.pi, [0, 0, 1]))
    sc.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [0, 1, 0]))
    sc.apply_transform(trimesh.transformations.rotation_matrix(np.radians(56), [0, 0, 1]))
    sc.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))
    sc.export(path)


def render(path_glb, png, samples):
    env = dict(os.environ, LC_ALL="C", LANG="C")
    r = subprocess.run([BLENDER, "-b", "-P", str(SCRIPT), "--", str(path_glb), str(png), "--samples", str(samples),
                        "--res", "1400", "--bg"], env=env, capture_output=True, text=True)
    if not Path(png).exists() or r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
        raise SystemExit("blender failed")


def main():
    out = HERE / "out"
    items = scene_items(0.5)
    glb(items, out / "wobble2.glb")
    glb(cutaway(items), out / "wobble2_cut.glb")
    n = 32 if "--quick" in sys.argv else 128
    render(out / "wobble2.glb", out / "render_full.png", n)
    render(out / "wobble2_cut.glb", out / "render_cut.png", n)
    core = [it for it in items if it[0] in ("ring", "rollers", "input_carrier", "output_carrier")]
    glb(core, out / "wobble2_core.glb")
    render(out / "wobble2_core.glb", out / "render_core.png", n)
    print("wrote out/render_full.png, out/render_cut.png, out/render_core.png, out/wobble2.glb")


if __name__ == "__main__":
    main()
