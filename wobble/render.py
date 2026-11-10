"""
Renders of the twin wobble reducer from cad.py's meshes (product-design skill pipeline:
assemble in world coords -> PBR GLB -> headless Blender).

    python wobble/render.py      (or: make wobble-render)   -> out/render_full.png, out/render_cut.png

The cutaway removes the half of the housings and output ring facing the camera (and leaves out
the clevis and motor) so the blocks, rollers and grooves show; moving parts stay whole.
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
STEEL = ((0.60, 0.62, 0.65, 1), 1.0, 0.25)
DARK = ((0.12, 0.12, 0.13, 1), 0.5, 0.45)
ORANGE = ((0.93, 0.45, 0.12, 1), 0.0, 0.45)
BLUE = ((0.20, 0.40, 0.70, 1), 0.3, 0.35)
BRG = ((0.82, 0.84, 0.87, 1), 1.0, 0.2)
CUT = ("housing_a", "housing_b", "output_ring")


def scene_items(phi=0.0):
    """(name, mesh in world, material) at input angle phi"""
    m = {p.stem: trimesh.load(p) for p in CAD.glob("*.stl")}
    items = [("housing_a", m["housing_a"], ALU), ("housing_b", m["housing_b"], ALU),
             ("output_ring", m["output_ring"], BLUE), ("clevis", m["clevis"], DARK), ("motor", m["motor"], DARK),
             ("input_shaft", m["input_shaft"], STEEL)]
    out = [(n, C.transformed(mesh, C.pose(n, phi)), mat) for n, mesh, mat in items]
    for k in (0, 1):
        out.append((f"block_{k}", C.transformed(m["block_canonical"], C.pose("block", phi, k)), ORANGE))
        out.append((f"wobble_bearing_{k}", C.transformed(m["wobble_bearing"], C.pose("wobble_bearing", phi, k)), BRG))
    rollers = [C.transformed(m["roller_canonical"], C.pose("roller", phi, k, u)) for _r, k, _j, u, _s in C.roller_axes()]
    out.append(("rollers", trimesh.util.concatenate(rollers), STEEL))
    return out


def cutaway(items):
    box = trimesh.creation.box(extents=(400, 400, 400))
    box.apply_translation((0, 200, 0))                       # the +y half (lug side, toward the camera)
    out = []
    for n, mesh, mat in items:
        if n in ("clevis", "motor"):
            continue
        if n in CUT:
            mesh = C.boolean(mesh, box, "diff")
        out.append((n, mesh, mat))
    return out


def glb(items, path):
    sc = trimesh.Scene()
    for n, mesh, (col, metal, rough) in items:
        mesh = mesh.copy()
        mesh.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=n, baseColorFactor=col, metallicFactor=metal, roughnessFactor=rough, doubleSided=True))
        sc.add_geometry(mesh, geom_name=n, node_name=n)
    # lug toward the camera (it looks from -y), axis horizontal (z -> x), turned a little, Z-up -> Y-up
    sc.apply_transform(trimesh.transformations.rotation_matrix(np.pi, [0, 0, 1]))
    sc.apply_transform(trimesh.transformations.rotation_matrix(np.pi / 2, [0, 1, 0]))
    sc.apply_transform(trimesh.transformations.rotation_matrix(np.radians(56), [0, 0, 1]))
    sc.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))
    sc.export(path)


def render(path_glb, png, samples=128):
    env = dict(os.environ, LC_ALL="C", LANG="C")
    r = subprocess.run([BLENDER, "-b", "-P", str(SCRIPT), "--", str(path_glb), str(png), "--samples", str(samples),
                        "--res", "1400", "--bg"], env=env, capture_output=True, text=True)
    if not Path(png).exists() or r.returncode:
        print(r.stdout[-2000:], r.stderr[-2000:])
        raise SystemExit("blender failed")


def main():
    out = HERE / "out"
    items = scene_items(0.6)
    glb(items, out / "wobble.glb")
    glb(cutaway(items), out / "wobble_cut.glb")
    quick = "--quick" in sys.argv
    render(out / "wobble.glb", out / "render_full.png", 32 if quick else 128)
    render(out / "wobble_cut.glb", out / "render_cut.png", 32 if quick else 128)
    print("wrote out/render_full.png, out/render_cut.png, out/wobble.glb")


if __name__ == "__main__":
    main()
