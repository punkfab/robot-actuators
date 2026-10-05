"""
Presentation renders of the shallow-slot motor: cad.py's assembly -> GLB with PBR
materials -> Blender (Cycles), using the product-design skill's blender_render.py.

Scenes
  cutaway   the whole motor with a window cut through the housing, a narrower one
            through the rotor, and coolant shown in the jacket channel
  stator    the wound stator on its own: core, 24 coils, manifold

Viz only: materials, cuts, tilts and the coolant body live here, not in cad.py.

    python qdd-liquid/render.py        (or: make qdd-liquid-render)
"""

import os
import subprocess
import sys
from math import cos, degrees, radians, sin
from pathlib import Path

import numpy as np
import trimesh
from build123d import Compound, Polygon, Pos, Rot, export_stl, extrude
from trimesh.visual.material import PBRMaterial

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import cad  # noqa: E402
from cad import D  # noqa: E402

OUT = HERE / "out" / "render"
BLENDER = os.environ.get("BLENDER", "/opt/blender-5.0.1-linux-x64/blender")
SCRIPT = Path.home() / ".claude/skills/product-design/scripts/blender_render.py"
CAMERA_AZ = -56.0          # where blender_render.py's camera sits, degrees in the XY plane

# name -> (rgba, metallic, roughness)
MATERIALS = {
    "stator_core": ((0.20, 0.22, 0.25, 1), 0.9, 0.45),
    "coil": ((0.85, 0.45, 0.22, 1), 1.0, 0.28),
    "rotor_yoke": ((0.16, 0.17, 0.19, 1), 0.9, 0.4),
    "magnet": ((0.78, 0.80, 0.82, 1), 1.0, 0.18),
    "rotor_carrier": ((0.62, 0.64, 0.67, 1), 1.0, 0.35),
    "jacket_inner": ((0.70, 0.72, 0.75, 1), 1.0, 0.35),
    "jacket_outer": ((0.10, 0.11, 0.13, 1), 0.6, 0.45),
    "coolant": ((0.10, 0.55, 0.85, 0.55), 0.0, 0.08),
    "manifold": ((0.80, 0.72, 0.55, 1), 0.0, 0.6),
    "drive_plate": ((0.12, 0.40, 0.28, 1), 0.3, 0.5),
    "cap_rear": ((0.10, 0.11, 0.13, 1), 0.6, 0.45),
    "cap_front": ((0.10, 0.11, 0.13, 1), 0.6, 0.45),
    "bearing": ((0.85, 0.86, 0.88, 1), 1.0, 0.15),
    "shaft": ((0.55, 0.57, 0.60, 1), 1.0, 0.3),
}


def wedge(a0, a1, r=95.0, z0=-15.0, h=110.0):
    """Prism over the sector a0..a1 (degrees), for cutting a window."""
    pts = [(0, 0)] + [(r * cos(radians(a)), r * sin(radians(a))) for a in np.linspace(a0, a1, 24)]
    return Pos(0, 0, z0) * extrude(Polygon(*pts, align=None), amount=h)


def minus(shape, tool):
    """shape - tool as a flat list of solids (a boolean can return one shape, several, or none)."""
    r = shape - tool
    parts = list(r) if isinstance(r, (list, tuple)) else [r]
    return [x for part in parts for x in (part.solids() if hasattr(part, "solids") else [part])]


def cut(solids, a0, a1):
    w = wedge(a0, a1)
    return [x for s in solids for x in minus(s, w)]


def cut_coils(coils, a0, a1):
    """Drop coils wholly inside the window, cut the ones straddling its edges."""
    w, out = wedge(a0, a1), []
    half = degrees(cad.SLOT_PITCH) / 2 + 1.0
    for k, c in enumerate(coils):
        a = degrees((k + 0.5) * cad.SLOT_PITCH)
        if a0 + half < a < a1 - half:
            continue
        if a < a0 - half or a > a1 + half:
            out.append(c)
        else:
            out += [x for s in c.solids() for x in minus(s, w)]
    return out


def coolant():
    """The oil in the jacket's annular channel."""
    h = cad.LENGTH - 2 * cad.CAP - 2 * cad.LAND
    return Pos(0, 0, cad.STATIONS["jacket"] + cad.LAND) * cad.ring(cad.R_STATOR + cad.J_WALL, cad.R_OUT - cad.J_WALL, h)


def scene_cutaway():
    p = cad.assembly()
    window, rotor_window = (-25, 125), (85, 125)
    s = {}
    for name in ("stator_core", "jacket_inner", "jacket_outer", "manifold", "drive_plate", "cap_rear", "cap_front"):
        s[name] = cut(p[name], *window)
    s["coolant"] = cut([coolant()], *window)
    s["coil"] = cut_coils(p["coil"], 60, 125)
    for name in ("rotor_yoke", "rotor_carrier", "bearing"):
        s[name] = cut(p[name], *rotor_window)
    s["magnet"] = [m for k, m in enumerate(p["magnet"])
                   if not rotor_window[0] < degrees(k * cad.POLE_PITCH) < rotor_window[1]]
    s["shaft"] = p["shaft"]
    spin = Rot(0, 0, CAMERA_AZ - 50.0)                     # window centre toward the camera
    return {k: [spin * x for x in v] for k, v in s.items()}


def scene_stator():
    p = cad.assembly()
    s = {k: p[k] for k in ("stator_core", "coil", "manifold")}
    # tip the front face toward the camera
    tilt = Rot(0, 0, CAMERA_AZ) * Rot(0, 38, 0) * Rot(0, 0, -CAMERA_AZ)
    z = cad.STATIONS["stack"] + cad.STACK / 2
    return {k: [tilt * (Pos(0, 0, -z) * x) for x in v] for k, v in s.items()}


def to_glb(name, parts):
    OUT.mkdir(parents=True, exist_ok=True)
    scene = trimesh.Scene()
    for part, solids in parts.items():
        stl = OUT / f"_{name}_{part}.stl"
        export_stl(Compound(children=list(solids)), str(stl), tolerance=0.02, angular_tolerance=0.12)
        mesh = trimesh.load(stl)
        rgba, metal, rough = MATERIALS[part]
        mesh.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=part, baseColorFactor=rgba, metallicFactor=metal, roughnessFactor=rough,
            alphaMode="BLEND" if rgba[3] < 1 else "OPAQUE", doubleSided=True))
        scene.add_geometry(mesh, geom_name=part)
        stl.unlink()
    scene.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))   # Z-up -> Y-up
    glb = OUT / f"{name}.glb"
    scene.export(glb)
    return glb


def blender(glb, png, samples=200, res=1400):
    if png.exists():
        png.unlink()
    env = dict(os.environ, LC_ALL="C", LANG="C")
    subprocess.run([BLENDER, "-b", "-P", str(SCRIPT), "--", str(glb), str(png), "--samples", str(samples),
                    "--res", str(res), "--bg"], env=env, check=True, capture_output=True)
    if not png.exists():
        raise RuntimeError(f"Blender wrote nothing for {png}")
    print(f"wrote {png}")


def main():
    for name, build in (("cutaway", scene_cutaway), ("stator", scene_stator)):
        glb = to_glb(name, build())
        print(f"{glb.name}: {glb.stat().st_size / 1e6:.1f} MB")
        blender(glb, OUT / f"{name}.png")


if __name__ == "__main__":
    main()
