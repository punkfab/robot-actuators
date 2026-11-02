"""
Meshes and motion data for the interactive assembly page (web/index.html): the same parts,
exploded layout and running kinematics as assemble_anim.py, but posed live in the browser.

    python rv/web_export.py          (or: make rv-web)  -> web/parts.gltf.json (glTF with embedded buffers), web/meta.json

then serve the folder (python -m http.server -d rv/web) and open http://localhost:8000
"""

import functools
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
WEB = HERE / "web"
sys.path.insert(0, str(HERE))
import assemble_anim  # noqa: E402
import explode  # noqa: E402
import printed  # noqa: E402

printed.build = functools.lru_cache(printed.build)          # objects() and groups() both build the CAD


def to_gltf(rows, path):
    """one node per part, CAD frame (mm, axis z); coarser tessellation than the render, for the web"""
    import trimesh
    from build123d import export_stl
    from trimesh.visual.material import PBRMaterial
    scene = trimesh.Scene()
    tmp = WEB / "_tmp.stl"
    for name, g, role, k, fade, col, solid in rows:
        export_stl(solid, str(tmp), tolerance=0.05, angular_tolerance=0.3)
        m = trimesh.load(tmp)
        m.merge_vertices()
        metal = 1.0 if col in (explode.BRG, explode.STEEL) else (0.6 if col == explode.BLACK else 0.0)
        m.visual = trimesh.visual.TextureVisuals(material=PBRMaterial(
            name=name, baseColorFactor=col, metallicFactor=metal, roughnessFactor=0.3 if metal else 0.55, doubleSided=True))
        scene.add_geometry(m, geom_name=name, node_name=name)
    tmp.unlink()
    from trimesh.exchange.gltf import export_gltf
    files = export_gltf(scene, embed_buffers=True)                 # one JSON file, buffers inline
    path.write_bytes(files["model.gltf"])


def main():
    WEB.mkdir(exist_ok=True)
    rows, order = assemble_anim.objects()
    dz = assemble_anim.explode_offsets(rows, order)
    to_gltf(rows, WEB / "parts.gltf.json")
    desc = {key: dict(text=text, qty=qty, kind=kind) for key, text, qty, kind, _c, _s in explode.groups()}
    P = printed.P
    meta = dict(
        ratio=P.ratio, e=P.e, r_crank=P.r_crank, crank0=printed.CRANK0, n_crank=P.n_crank,
        z_sun=P.z_sun, z_planet=P.z_planet, pins=P.pins, lobes=P.lobes,
        groups=[dict(key=key, name=assemble_anim.GROUP_NAMES[g], dz=dz[g], **desc.get(key, {}))
                for g, key in enumerate(order)],
        objects=[dict(name=r[0], group=r[1], role=r[2], k=r[3], fade=r[4]) for r in rows])
    (WEB / "meta.json").write_text(json.dumps(meta, indent=1))
    size = (WEB / "parts.gltf.json").stat().st_size / 1e6
    print(f"wrote {WEB}/parts.gltf.json ({len(rows)} parts, {size:.1f} MB) and meta.json ({len(order)} assembly steps)")


if __name__ == "__main__":
    main()
