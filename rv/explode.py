"""
Exploded assembly drawing of the printable RV (printed.py): every part group pulled apart
along the axis in assembly order (motor end first), rendered in Blender, then numbered
balloons with leader lines and a parts list drawn on with PIL.

    python rv/explode.py          (or: make rv-explode)  -> out/explode_printed.png
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out" / "explode"
sys.path.insert(0, str(HERE))

BLENDER = os.environ.get("BLENDER", "/opt/blender-5.0.1-linux-x64/blender")
GAP = 26.0                         # mm between exploded groups along the axis
W, H, PANEL = 2800, 1250, 1000     # render size and the parts-list panel to its right
FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"

ORANGE = (0.95, 0.45, 0.10, 1)
BLUE = (0.18, 0.38, 0.68, 1)
WHITE = (0.88, 0.88, 0.85, 1)
STEEL = (0.62, 0.64, 0.67, 1)
BRG = (0.80, 0.82, 0.85, 1)
BLACK = (0.10, 0.10, 0.11, 1)


def groups():
    """(key, part-list text, qty, bought/printed, colour, solids) in assembly order"""
    import printed
    p, _ = printed.build()
    eb = p["ecc_bearing"]                       # per crank: disc 1 (x2), disc 2 (x2)
    jb = p["journal_bearing"]                   # per crank: rear, front
    n = len(p["crank_a"])
    return [
        ("motor", "NEMA 17 stepper, 5 mm D shaft", 1, "bought", BLACK, p["motor"]),
        ("input_cover", "input cover / motor mount", 1, "printed", WHITE, p["input_cover"]),
        ("sun", "sun gear, 12 teeth m1, on the motor shaft", 1, "printed", BLUE, p["sun"]),
        ("crank_a", "crankshaft piece A: planet (32 t) + rear journal + hex", n, "printed", ORANGE, p["crank_a"]),
        ("main_rear", "main bearing 6814-2RS, 70x90x10", 1, "bought", BRG, [p["main_bearing"][0]]),
        ("journal_rear", "journal bearing 6800-2RS, 10x19x5", n, "bought", BRG, [jb[2 * k] for k in range(n)]),
        ("carrier_rear", "carrier plate, rear", 1, "printed", BLUE, p["carrier_rear"]),
        ("ecc1", "eccentric bearings 6802-2RS, 15x24x5 (2 per crank)", 2 * n, "bought", BRG,
         [eb[4 * k + j] for k in range(n) for j in (0, 1)]),
        ("disc1", "cycloid disc (29 lobes)", 1, "printed", ORANGE, [p["disc"][0]]),
        ("ecc2", "eccentric bearings 6802-2RS, 15x24x5 (2 per crank)", 2 * n, "bought", BRG,
         [eb[4 * k + j] for k in range(n) for j in (2, 3)]),
        ("disc2", "cycloid disc, the same part flipped over", 1, "printed", ORANGE, [p["disc"][1]]),
        ("crank_b", "crankshaft piece B: both eccentrics + front journal", n, "printed", ORANGE, p["crank_b"]),
        ("crank_core", "crankshaft core, DIN 6325 dowel 5x45", n, "bought", STEEL, p["crank_core"]),
        ("pillar", "carrier pillar, DIN 6325 dowel 5x45", n, "bought", STEEL, p["pillar"]),
        ("pin", "ring pin, DIN 6325 dowel 5x24", len(p["pin"]), "bought", STEEL, p["pin"]),
        ("housing", "housing / pin case", 1, "printed", WHITE, p["housing"]),
        ("journal_front", "journal bearing 6800-2RS, 10x19x5", n, "bought", BRG, [jb[2 * k + 1] for k in range(n)]),
        ("main_front", "main bearing 6814-2RS, 70x90x10", 1, "bought", BRG, [p["main_bearing"][1]]),
        ("carrier_front", "carrier plate, front (output flange)", 1, "printed", BLUE, p["carrier_front"]),
    ]


def exploded(gs):
    """stack the groups along z: each starts GAP after the previous one ends"""
    from build123d import Compound, Pos
    out, z_next = [], None
    for key, text, qty, kind, col, solids in gs:
        bb = Compound(children=solids).bounding_box()
        dz = 0.0 if z_next is None else z_next - bb.min.Z
        moved = [Pos(0, 0, dz) * s for s in solids]
        out.append((key, text, qty, kind, col, moved))
        z_next = bb.max.Z + dz + GAP
    return out


def to_glb(gs, path):
    import trimesh
    from build123d import Compound, Rot, export_stl
    from trimesh.visual.material import PBRMaterial
    scene = trimesh.Scene()
    lay = Rot(0, 90, 0)                                  # axis z -> x
    for key, text, qty, kind, col, solids in gs:
        stl = OUT / f"_{key}.stl"
        export_stl(Compound(children=[lay * s for s in solids]), str(stl), tolerance=0.02, angular_tolerance=0.15)
        m = trimesh.load(stl)
        metal = 1.0 if kind == "bought" and col != BLACK else (0.6 if col == BLACK else 0.0)
        m.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=key, baseColorFactor=col, metallicFactor=metal, roughnessFactor=0.3 if metal else 0.6, doubleSided=True))
        scene.add_geometry(m, geom_name=key, node_name=key)
        stl.unlink()
    scene.apply_transform(trimesh.transformations.rotation_matrix(-np.pi / 2, [1, 0, 0]))
    scene.export(path)


def annotate(gs, png, labels, out_png, title):
    from PIL import Image, ImageDraw, ImageFont
    img = Image.open(png).convert("RGB")
    canvas = Image.new("RGB", (W + PANEL, H), (246, 246, 247))
    canvas.paste(img, (0, 0))
    dr = ImageDraw.Draw(canvas)
    f_b = ImageFont.truetype(FONT, 26)
    f = ImageFont.truetype(FONT.replace("Sans.ttf", "Sans.ttf"), 22)
    f_t = ImageFont.truetype(FONT.replace("Sans.ttf", "Sans-Bold.ttf"), 30)
    f_h = ImageFont.truetype(FONT.replace("Sans.ttf", "Sans-Bold.ttf"), 22)
    # balloons: alternate rows above and below, spread so they don't overlap
    items = []
    for i, (key, *_rest) in enumerate(gs):
        lab = labels[key]
        above = i % 2 == 0
        anchor = lab["center"]                          # a vertex on the part, nearest its projected centre
        items.append(dict(n=i + 1, above=above, ax=anchor[0], ay=anchor[1]))
    r = 26
    for above in (True, False):
        row = sorted([it for it in items if it["above"] == above], key=lambda it: it["ax"])
        xs = [it["ax"] for it in row]
        for k in range(1, len(xs)):                      # push right to keep 2.6 r apart
            xs[k] = max(xs[k], xs[k - 1] + 2.6 * r)
        shift = max(0.0, xs[-1] - (W - 2 * r)) if xs else 0
        xs = [x - shift for x in xs]
        for k in range(len(xs) - 2, -1, -1):             # then pull back left if pushed off the edge
            xs[k] = min(xs[k], xs[k + 1] - 2.6 * r)
        for it, x in zip(row, xs):
            y = 95 if above else H - 70
            it["bx"], it["by"] = max(x, r + 4), y
    for it in items:
        bx, by = it["bx"], it["by"]
        ey = by + (r if it["above"] else -r)
        dr.line([(bx, ey), (it["ax"], it["ay"])], fill=(90, 90, 95), width=2)
        dr.ellipse([it["ax"] - 4, it["ay"] - 4, it["ax"] + 4, it["ay"] + 4], fill=(90, 90, 95))
        dr.ellipse([bx - r, by - r, bx + r, by + r], fill=(255, 255, 255), outline=(40, 40, 45), width=3)
        t = str(it["n"])
        w = dr.textlength(t, font=f_b)
        dr.text((bx - w / 2, by - 15), t, fill=(20, 20, 25), font=f_b)
    dr.text((30, 20), title, fill=(20, 20, 25), font=f_t)
    # parts list
    x0, y0 = W + 30, 40
    dr.text((x0, y0), "Parts list", fill=(20, 20, 25), font=f_t)
    y = y0 + 56
    dr.text((x0, y), "#", fill=(60, 60, 65), font=f_h)
    dr.text((x0 + 50, y), "part", fill=(60, 60, 65), font=f_h)
    dr.text((x0 + 760, y), "qty", fill=(60, 60, 65), font=f_h)
    dr.text((x0 + 830, y), "", fill=(60, 60, 65), font=f_h)
    y += 36
    for i, (key, text, qty, kind, col, _s) in enumerate(gs):
        if i % 2 == 0:
            dr.rectangle([x0 - 10, y - 4, W + PANEL - 20, y + 48], fill=(236, 236, 239))
        sw = tuple(int(c * 255) for c in col[:3])
        dr.rectangle([x0 + 870, y + 6, x0 + 900, y + 30], fill=sw, outline=(90, 90, 95))
        dr.text((x0, y + 8), str(i + 1), fill=(20, 20, 25), font=f)
        dr.text((x0 + 50, y), text, fill=(20, 20, 25), font=f)
        dr.text((x0 + 50, y + 24), kind, fill=(120, 120, 125) if kind == "printed" else (40, 90, 160), font=ImageFont.truetype(FONT, 17))
        dr.text((x0 + 770, y + 8), str(qty), fill=(20, 20, 25), font=f)
        y += 56
    canvas.save(out_png)


def main():
    import printed
    OUT.mkdir(parents=True, exist_ok=True)
    gs = exploded(groups())
    glb = OUT / "printed.glb"
    to_glb(gs, glb)
    png, js = OUT / "_render.png", OUT / "labels.json"
    env = dict(os.environ, LC_ALL="C", LANG="C")
    r = subprocess.run([BLENDER, "-b", "-P", str(HERE / "explode_blender.py"), "--", str(glb), str(png), str(js),
                        str(W), str(H), "96"], env=env, capture_output=True, text=True)
    if not png.exists() or r.returncode:
        print(r.stdout[-3000:], r.stderr[-3000:])
        raise SystemExit("blender failed")
    labels = json.loads(js.read_text())
    title = (f"Printable RV reducer, {printed.P.ratio:.0f}:1 - exploded view "
             f"({sum(g[2] for g in gs)} parts: {sum(g[2] for g in gs if g[3] == 'printed')} printed)")
    out_png = HERE / "out" / "explode_printed.png"
    annotate(gs, png, labels, out_png, title)
    print(f"wrote {out_png}")


if __name__ == "__main__":
    main()
