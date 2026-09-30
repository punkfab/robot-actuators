"""
RobStride CAD for assemblies: the vendor STEP, in one frame for every model.

RobStride publishes a STEP per model (github.com/RobStride/Product_Information,
"Product Literature/"). No licence is stated, so they're fetched on demand into
robstride/cad/ (git-ignored), not vendored.

Each vendor file has its own origin and orientation. actuator(kind) returns it
normalized to the frame every design here uses:

    +Z   the joint axis, pointing out of the output flange
    z=0  the output flange face (bolt your link here; the housing is at z < 0)
    x,y  centred on the axis

The axis is the bbox direction whose two perpendicular extents match (round
housing). The output end is the end the output flange's tapped holes open from:
M3–M5 holes on the booklet's output bolt circle (Ø24–30), which is 2–3× smaller than
the housing's. `python cad.py` checks every model against the booklet's OD×H and
the expected end, and writes a line-up STEP.

envelope(kind) is the no-network stand-in: a cylinder of the booklet's OD×H.
"""

import sys
import urllib.request
from pathlib import Path

import numpy as np
from build123d import Cylinder, Location, Pos, Rot, export_brep, import_brep, import_step

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from specs import SPECS  # noqa: E402

CACHE = HERE / "cad"
BASE = "https://raw.githubusercontent.com/RobStride/Product_Information/main/Product%20Literature/"
FILES = {"RS00": "RS00/RS00.step", "RS01": "RS01/RS01.stp", "RS02": "RS02/rs02.stp",
         "RS03": "RS03/RS03.stp", "RS04": "RS04/RS04.stp", "RS05": "RS05/RS05.STEP",
         "RS06": "RS06/RS06-new.step", "RS10P": "10P/RS10.stp"}


def fetch(kind) -> Path:
    CACHE.mkdir(exist_ok=True)
    path = CACHE / f"{kind}{Path(FILES[kind]).suffix.lower()}"
    if not path.exists():
        urllib.request.urlretrieve(BASE + FILES[kind], path)
    return path


# Output bolt circle (PCD, mm) from the booklet drawings [SPEC pp. 2–34]. The
# output end is the end these tapped holes open from; the housing circles are 2–3×
# larger. RS01 and RS05 [INFERRED]: RS01 from RS02's drawing, RS05's callout ambiguous.
OUTPUT_PCD = {"RS00": 27, "RS01": 24, "RS02": 27, "RS03": 30.36, "RS04": 30.36,
              "RS05": 24, "RS06": 24, "RS10P": 27}
# What output_end() finds (checked by `python cad.py`): +1 = vendor's +axis end.
EXPECTED_END = {"RS00": -1, "RS01": -1, "RS02": 1, "RS03": 1, "RS04": 1, "RS05": 1,
                "RS06": -1, "RS10P": -1}


def _cyl_radius(f):
    from OCP.BRepAdaptor import BRepAdaptor_Surface
    return BRepAdaptor_Surface(f.wrapped).Cylinder().Radius()


def output_end(shape, pcd):
    """+1/−1: which z end the output bolt circle opens from (shape axis-aligned on z,
    centred). Uses M3–M5-sized cylindrical faces whose centres sit on the PCD."""
    from build123d import GeomType
    zs = []
    for f in shape.faces():
        if f.geom_type != GeomType.CYLINDER or not (1.2 <= _cyl_radius(f) <= 3.2):
            continue
        p = f.center()
        if abs(np.hypot(p.X, p.Y) - pcd / 2) < 0.8:
            zs.append(p.Z)
    if not zs:
        raise ValueError("no holes on the output bolt circle")
    return 1 if np.mean(zs) > 0 else -1


def _normalize(shape, pcd):
    """Rotate/translate the vendor shape into the frame above. Returns (shape, info)."""
    bb = shape.bounding_box()
    size = np.array([bb.size.X, bb.size.Y, bb.size.Z])
    # axis: the one whose two perpendicular extents are closest to equal
    pairs = [abs(size[1] - size[2]), abs(size[0] - size[2]), abs(size[0] - size[1])]
    ax = int(np.argmin(pairs))
    s = {0: Rot(0, 90, 0), 1: Rot(-90, 0, 0), 2: Rot(0, 0, 0)}[ax] * shape
    c = s.bounding_box().center()
    s = Pos(-c.X, -c.Y, -c.Z) * s
    end = output_end(s, pcd)
    if end < 0:
        s = Rot(180, 0, 0) * s
    s = Pos(0, 0, -s.bounding_box().max.Z) * s
    bb = s.bounding_box()
    return s, dict(od=max(bb.size.X, bb.size.Y), h=bb.size.Z, end=end)


_cache = {}


def raw(kind):
    """Vendor shape as published (STEP import is slow: cached as BREP next to it)."""
    brep = CACHE / f"{kind}.brep"
    if not brep.exists():
        export_brep(import_step(str(fetch(kind))), str(brep))
    return import_brep(str(brep))


def actuator(kind, loc: Location | None = None):
    """Vendor STEP of `kind`, normalized (see module doc); optionally placed at loc.
    The normalized shape is cached as cad/<kind>_norm.brep."""
    if kind not in _cache:
        norm = CACHE / f"{kind}_norm.brep"
        if norm.exists():
            s = import_brep(str(norm))
            bb = s.bounding_box()
            _cache[kind] = (s, dict(od=max(bb.size.X, bb.size.Y), h=bb.size.Z,
                                    end=EXPECTED_END[kind]))
        else:
            _cache[kind] = _normalize(raw(kind), OUTPUT_PCD[kind])
            export_brep(_cache[kind][0], str(norm))
    s = _cache[kind][0]
    return loc * s if loc is not None else s


def envelope(kind, loc: Location | None = None):
    """Booklet OD×H cylinder in the same frame (housing below z = 0)."""
    sp = SPECS[kind]
    c = Pos(0, 0, -sp.length * 500) * Cylinder(sp.od * 500, sp.length * 1000)
    return loc * c if loc is not None else c


def main():
    from build123d import Compound, export_step
    print("model   STEP OD×H (mm)      booklet OD×H    output end (vendor axis)")
    parts, x = [], 0.0
    for kind in FILES:
        s = actuator(kind)
        info = _cache[kind][1]
        sp = SPECS[kind]
        print(f"{kind:6s}  {info['od']:6.1f} × {info['h']:5.1f}      "
              f"{sp.od*1e3:5.1f} × {sp.length*1e3:5.1f}    "
              f"{'+' if info['end'] > 0 else '−'}z  "
              f"{'ok' if info['end'] == EXPECTED_END[kind] else 'CHANGED'}")
        x += info["od"] / 2
        parts.append(Pos(x, 0, 0) * s)
        x += info["od"] / 2 + 20
    out = HERE / "out"
    out.mkdir(exist_ok=True)
    export_step(Compound(children=parts), str(out / "lineup.step"))
    print(f"wrote {out/'lineup.step'} (all models, output flange up at z = 0)")


if __name__ == "__main__":
    main()
