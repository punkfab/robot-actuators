"""
Presentation render of the RV actuator: cad.py's assembly -> GLB with PBR materials ->
Blender (Cycles), using the product-design skill's blender_render.py.

  section   half-section through the axis, lying along x: motor, planetary stage,
            crankshafts with their needle rollers, both cycloid discs on the pins, carrier

Viz only: materials and the cut live here, not in cad.py.

    python rv/render.py        (or: make rv-render)
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

STEEL = ((0.66, 0.68, 0.71, 1), 1.0, 0.28)
DARK = ((0.12, 0.13, 0.15, 1), 0.6, 0.45)
BRG = ((0.85, 0.86, 0.88, 1), 1.0, 0.15)
MATERIALS = {
    "housing": ((0.72, 0.74, 0.77, 1), 1.0, 0.38), "pin": ((0.85, 0.86, 0.88, 1), 1.0, 0.12),
    "disc": ((0.80, 0.62, 0.30, 1), 1.0, 0.30), "crankshaft": ((0.55, 0.57, 0.60, 1), 1.0, 0.25),
    "roller": ((0.90, 0.90, 0.92, 1), 1.0, 0.10), "ecc_ring": ((0.80, 0.82, 0.85, 1), 1.0, 0.15), "carrier_front": ((0.35, 0.45, 0.60, 1), 0.9, 0.35),
    "carrier_rear": ((0.35, 0.45, 0.60, 1), 0.9, 0.35), "journal_bearing": BRG, "main_bearing": BRG,
    "sun": ((0.75, 0.60, 0.35, 1), 1.0, 0.3), "planet": ((0.75, 0.60, 0.35, 1), 1.0, 0.3), "input_shaft": STEEL,
    "input_bearing": BRG, "gear_cover": DARK, "motor_can": DARK, "stator": ((0.25, 0.27, 0.30, 1), 0.9, 0.45),
    "end_windings": ((0.85, 0.45, 0.22, 1), 1.0, 0.28), "rotor": ((0.20, 0.20, 0.22, 1), 0.9, 0.4),
}


def solids_of(r):
    parts = list(r) if isinstance(r, (list, tuple)) else [r]
    return [x for p in parts for x in (p.solids() if hasattr(p, "solids") else [p])]


KEEP = Pos(0, 0, -50) * Box(500, 250, 500, align=(Align.CENTER, Align.MIN, Align.MIN))     # y >= 0


def cut(solids):
    out = []
    for s in solids:
        bb = s.bounding_box()
        if bb.min.Y >= 0:
            out.append(s)
        elif bb.max.Y > 0:
            out += solids_of(s & KEEP)
    return out


LAY = Rot(0, 90, 0)          # axis z -> x

PLA_A = ((0.95, 0.45, 0.10, 1), 0.0, 0.55)          # printed parts: matte filament
PLA_B = ((0.15, 0.35, 0.65, 1), 0.0, 0.55)
PLA_C = ((0.85, 0.85, 0.82, 1), 0.0, 0.6)
MATERIALS.update({"ecc_bearing": BRG, "pillar": STEEL, "crank_core": STEEL, "crank_a": PLA_A, "crank_b": PLA_A,
                  "input_cover": PLA_C, "motor": ((0.10, 0.10, 0.11, 1), 0.7, 0.4)})
PRINTED_MAT = {"housing": PLA_C, "disc": PLA_A, "carrier_rear": PLA_B, "carrier_front": PLA_B, "sun": PLA_B}


def scene_printed():
    import printed
    p, _ = printed.build()
    turn = Rot(0, 0, -printed.CRANK0)
    p = {k: [turn * x for x in v] for k, v in p.items()}
    s = {}
    for k, v in p.items():
        if k == "pin":
            s[k] = [x for x in v if x.bounding_box().center().Y > -0.5]
        else:
            s[k] = cut(v)
    return {k: [LAY * x for x in v] for k, v in s.items() if v}


def scene_section():
    cad.phase_by_trial()
    turn = Rot(0, 0, -cad.CRANK0)                  # crank 1 onto +x, so the cut runs through it
    p = {k: [turn * x for x in v] for k, v in cad.assembly().items()}
    s = {}
    for k, v in p.items():
        if k in ("roller", "pin"):
            s[k] = [x for x in v if x.bounding_box().center().Y > -0.5]
        else:
            s[k] = cut(v)
    return {k: [LAY * x for x in v] for k, v in s.items() if v}


def to_glb(name, parts):
    OUT.mkdir(parents=True, exist_ok=True)
    scene = trimesh.Scene()
    for part, solids in parts.items():
        stl = OUT / f"_{name}_{part}.stl"
        export_stl(Compound(children=list(solids)), str(stl), tolerance=0.01, angular_tolerance=0.1)
        mesh = trimesh.load(stl)
        rgba, metal, rough = (PRINTED_MAT[part] if name == "printed" and part in PRINTED_MAT else MATERIALS[part])
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
    which = sys.argv[1:] or ["section", "printed"]
    for name, build in (("section", scene_section), ("printed", scene_printed)):
        if name not in which:
            continue
        glb = to_glb(name, build())
        print(f"{glb.name}: {glb.stat().st_size / 1e6:.1f} MB")
        blender(glb, OUT / f"{name}.png")


if __name__ == "__main__":
    main()
