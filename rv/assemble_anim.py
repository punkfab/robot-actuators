"""
Assembly-and-run animation of the printable RV (printed.py), rendered in Blender:

  act 1  the parts start in explode.py's exploded layout and slide home one group at a time,
         in assembly order, while the camera closes in
  act 2  housing, cover and main bearings fade to glass and the reducer runs: the input sun
         spins, the crankshafts turn and orbit with the carrier, the discs wobble on the pins,
         the output plate creeps round (the same kinematics as sim.py, checked there)

    python rv/assemble_anim.py          (or: make rv-assemble-video)  -> out/rv_assembly.mp4
    python rv/assemble_anim.py --quick  low resolution and samples, to check the motion
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "assemble"
sys.path.insert(0, str(HERE))
import explode  # noqa: E402

BLENDER = os.environ.get("BLENDER", "/opt/blender-5.0.1-linux-x64/blender")
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"


def objects():
    """every part as its own object: name, explode group, kinematic role, crank index, fade"""
    import printed
    p, _ = printed.build()
    n = len(p["crank_a"])
    eb, jb = p["ecc_bearing"], p["journal_bearing"]
    rows = []                        # (name, group, role, k, fade, colour, solid)

    def add(name, group, role, k, solid, colour, fade=False):
        rows.append((name, group, role, k, fade, colour, solid))

    order = ["motor", "input_cover", "sun", "crank_a", "main_rear", "journal_rear", "carrier_rear", "ecc1", "disc1", "ecc2",
             "disc2", "crank_b", "crank_core", "pillar", "pin", "housing", "journal_front", "main_front", "carrier_front"]
    gi = {g: i for i, g in enumerate(order)}
    add("motor", gi["motor"], "fixed", 0, p["motor"][0], explode.BLACK)
    add("input_cover", gi["input_cover"], "fixed", 0, p["input_cover"][0], explode.WHITE, fade=True)
    add("sun", gi["sun"], "input", 0, p["sun"][0], explode.BLUE)
    for k in range(n):
        add(f"crank_a_{k}", gi["crank_a"], "crank", k, p["crank_a"][k], explode.ORANGE)
        add(f"journal_rear_{k}", gi["journal_rear"], "carrier", 0, jb[2 * k], explode.BRG)
        add(f"journal_front_{k}", gi["journal_front"], "carrier", 0, jb[2 * k + 1], explode.BRG)
        add(f"crank_b_{k}", gi["crank_b"], "crank", k, p["crank_b"][k], explode.ORANGE)
        add(f"crank_core_{k}", gi["crank_core"], "crank", k, p["crank_core"][k], explode.STEEL)
        add(f"pillar_{k}", gi["pillar"], "carrier", 0, p["pillar"][k], explode.STEEL)
        for j in (0, 1):
            add(f"ecc1_{k}_{j}", gi["ecc1"], "disc1", k, eb[4 * k + j], explode.BRG)
            add(f"ecc2_{k}_{j}", gi["ecc2"], "disc2", k, eb[4 * k + 2 + j], explode.BRG)
    add("main_rear", gi["main_rear"], "fixed", 0, p["main_bearing"][0], explode.BRG, fade=True)
    add("main_front", gi["main_front"], "fixed", 0, p["main_bearing"][1], explode.BRG, fade=True)
    add("carrier_rear", gi["carrier_rear"], "carrier", 0, p["carrier_rear"][0], explode.BLUE)
    add("carrier_front", gi["carrier_front"], "carrier", 0, p["carrier_front"][0], explode.BLUE)
    add("disc1", gi["disc1"], "disc1", 0, p["disc"][0], explode.ORANGE)
    add("disc2", gi["disc2"], "disc2", 0, p["disc"][1], (0.85, 0.35, 0.08, 1))
    for i, pin in enumerate(p["pin"]):
        add(f"pin_{i}", gi["pin"], "fixed", 0, pin, explode.STEEL)
    add("housing", gi["housing"], "fixed", 0, p["housing"][0], explode.WHITE, fade=True)
    return rows, order


def explode_offsets(rows, order):
    """the same stacking as explode.py: each group starts GAP after the previous one ends"""
    from build123d import Compound
    dz, z_next = {}, None
    for g, _name in enumerate(order):
        solids = [r[6] for r in rows if r[1] == g]
        bb = Compound(children=solids).bounding_box()
        dz[g] = 0.0 if z_next is None else z_next - bb.min.Z
        z_next = bb.max.Z + dz[g] + explode.GAP
    return dz


def to_glb(rows, path):
    import trimesh
    from build123d import export_stl
    from trimesh.visual.material import PBRMaterial
    scene = trimesh.Scene()
    tmp = OUT / "_tmp.stl"
    for name, g, role, k, fade, col, solid in rows:
        export_stl(solid, str(tmp), tolerance=0.02, angular_tolerance=0.15)
        m = trimesh.load(tmp)
        metal = 1.0 if col in (explode.BRG, explode.STEEL) else (0.6 if col == explode.BLACK else 0.0)
        m.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=name, baseColorFactor=col, metallicFactor=metal, roughnessFactor=0.3 if metal else 0.55, doubleSided=True))
        scene.add_geometry(m, geom_name=name, node_name=name)      # CAD frame (mm, axis z); the Blender side maps it
    tmp.unlink()
    scene.export(path)


def captions(frames_dir, meta):
    from PIL import Image, ImageDraw, ImageFont
    f_b = ImageFont.truetype(FONT.replace("Sans.ttf", "Sans-Bold.ttf"), meta["font"])
    f = ImageFont.truetype(FONT, int(meta["font"] * 0.8))
    for png in sorted(frames_dir.glob("f*.png")):
        i = int(png.stem[1:])
        im = Image.open(png).convert("RGB")
        dr = ImageDraw.Draw(im)
        dr.text((30, 22), f"Printable RV reducer, {meta['ratio']:.0f}:1", fill=(25, 25, 30), font=f_b)
        if i < meta["act2"]:
            g = min(meta["n_groups"] - 1, max(0, (i - meta["t0"]) // meta["stagger"]))
            dr.text((30, 22 + meta["font"] * 1.4), f"assembly: {g + 1}/{meta['n_groups']}  {meta['names'][g]}",
                    fill=(70, 70, 75), font=f)
        else:
            s = meta["state"][i - meta["act2"]] if i - meta["act2"] < len(meta["state"]) else meta["state"][-1]
            dr.text((30, 22 + meta["font"] * 1.4), f"input {s[0]:6.2f} rev   crankshafts {s[1]:5.2f} rev   output {s[2]:5.1f} deg",
                    fill=(70, 70, 75), font=f)
        im.save(png)


def main():
    quick = "--quick" in sys.argv
    import printed
    OUT.mkdir(parents=True, exist_ok=True)
    rows, order = objects()
    dz = explode_offsets(rows, order)
    glb = OUT / "parts.glb"
    to_glb(rows, glb)
    P = printed.P
    fps = 30
    stagger, move, hold, run = 12, 26, 20, 240
    t0 = 15
    act2 = t0 + (len(order) - 1) * stagger + move + hold
    crank_revs = 3.0
    meta = dict(objects=[dict(name=r[0], group=r[1], role=r[2], k=r[3], fade=r[4]) for r in rows], dz=dz,
                e=P.e, r_crank=P.r_crank, crank0=printed.CRANK0, n_crank=P.n_crank, z_sun=P.z_sun, z_planet=P.z_planet,
                ratio=P.ratio, t0=t0, stagger=stagger, move=move, act2=act2, run=run, fps=fps, crank_revs=crank_revs,
                n_groups=len(order), quick=quick, res=(960, 540) if quick else (1920, 1080), samples=12 if quick else 48,
                font=22 if quick else 40,
                names=["NEMA 17 motor", "input cover", "sun gear", "crankshafts, piece A (planets)", "rear main bearing",
                       "rear journal bearings", "rear carrier plate", "eccentric bearings", "cycloid disc 1",
                       "eccentric bearings", "cycloid disc 2 (flipped)", "crankshafts, piece B (eccentrics)",
                       "steel crankshaft cores", "steel carrier pillars", "30 steel ring pins", "housing",
                       "front journal bearings", "front main bearing", "output carrier plate"])
    # running state for the captions
    w_c = crank_revs * 2 * np.pi / run
    state = []
    for f in range(run + 1):
        th_c = w_c * f
        th_in = th_c / (P.z_sun / P.z_planet) / (1 - 1 / P.ratio)
        state.append((th_in / (2 * np.pi), th_c / (2 * np.pi), np.degrees(th_in / P.ratio)))
    meta["state"] = state
    js = OUT / "meta.json"
    js.write_text(json.dumps(meta))
    frames = OUT / "frames"
    frames.mkdir(exist_ok=True)
    for f in frames.glob("*.png"):
        f.unlink()
    env = dict(os.environ, LC_ALL="C", LANG="C")
    r = subprocess.run([BLENDER, "-b", "-P", str(HERE / "assemble_blender.py"), "--", str(glb), str(js), str(frames)],
                       env=env, capture_output=True, text=True)
    n = len(list(frames.glob("f*.png")))
    if r.returncode or n == 0:
        print(r.stdout[-4000:], r.stderr[-4000:])
        raise SystemExit("blender failed")
    captions(frames, meta)
    mp4 = HERE / "out" / ("rv_assembly_quick.mp4" if quick else "rv_assembly.mp4")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", str(frames / "f%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "18", str(mp4)], check=True)
    print(f"wrote {mp4} ({n} frames, {mp4.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
