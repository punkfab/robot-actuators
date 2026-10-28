"""
Kinematic MuJoCo playback of the RV actuator: every moving part of cad.py's assembly as a
mocap body, posed each frame from the reducer's exact kinematics. It shows the motion, not
the dynamics: no forces, contact or friction are solved here.

  input (motor rotor, shaft, sun)   theta_in
  carrier (output)                  theta_o = theta_in / i,  i = 1 + (z_planet / z_sun) z_pins
  crankshafts and planets           relative to the carrier: theta_c = -(z_sun / z_planet)(theta_in - theta_o);
                                    each also orbits with the carrier
  cycloid discs                     turn with the carrier (theta_o); centres orbit at e, following
                                    the eccentrics (disc 2 opposite)
  needle rollers                    orbit their eccentric at about half the crank's relative speed
  housing, pins, motor stator       fixed

Before rendering, check() poses the discs at many input angles and measures, in 2-D, each
disc's clearance to all 40 pins: if any rate or sign were wrong the discs would cut the pins.

    python rv/sim.py              # interactive viewer, slowed down
    python rv/sim.py --video      # out/rv_sim.mp4 and out/rv_sim.gif
    python rv/sim.py --check      # the meshing check only
"""

import os
import sys

_MODE = "video" if "--video" in sys.argv else "check" if "--check" in sys.argv else "interactive"
os.environ.setdefault("MUJOCO_GL", "glfw" if _MODE == "interactive" else "egl")

import json  # noqa: E402
import subprocess  # noqa: E402
from math import cos, radians, sin, tau  # noqa: E402
from pathlib import Path  # noqa: E402

import mujoco  # noqa: E402
import numpy as np  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
MESH = OUT / "sim_meshes"
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402

J = json.loads((OUT / "design.json").read_text())
RV = dz.RV(**J["rv"])
I = RV.ratio
E = RV.e
CRANK0 = 90.0                                   # as cad.py


def crank_xy(k):
    a = radians(CRANK0) + tau * k / RV.n_crank
    return np.array([RV.r_crank * cos(a), RV.r_crank * sin(a)])


def angles(theta_in):
    th_o = theta_in / I
    th_c = -(RV.z_sun / RV.z_planet) * (theta_in - th_o)
    return th_o, th_c


def rot(a):
    return np.array([[cos(a), -sin(a)], [sin(a), cos(a)]])


# ---- poses: (rotation angle about z, xy translation) applied to meshes in their build pose -----
def pose_spin(a):
    return a, np.zeros(2)


def pose_crank(k, th_o, th_c):
    p = crank_xy(k)
    return th_o + th_c, rot(th_o) @ p - rot(th_o + th_c) @ p


def disc_centre(th_o, th_c):
    phi = th_o + th_c                            # each eccentric's direction in the world
    return E * np.array([cos(phi), sin(phi)])


def pose_disc(side, th_o, th_c):
    """side +1: disc 1 (built at (+e, 0)); -1: disc 2 (built at (-e, 0))"""
    c = side * disc_centre(th_o, th_c)
    built = np.array([side * E, 0.0])
    return th_o, c - rot(th_o) @ built


def pose_rollers(k, side, th_o, th_c):
    q0 = crank_xy(k) + np.array([side * E, 0.0])
    q = rot(th_o) @ crank_xy(k) + side * rot(th_o + th_c) @ np.array([E, 0.0])
    a = th_o + 0.5 * th_c                        # cage: about half the eccentric's speed relative to the disc
    return a, q - rot(a) @ q0


# ---- the meshing check -------------------------------------------------------------------------
def check(n=240, revs_in=3.0):
    from shapely import affinity
    from shapely.geometry import Point, Polygon
    x, y = dz.disc_profile(RV, steps=4000, mod=0.03)
    disc = Polygon(np.c_[x, y])
    disc2 = affinity.scale(disc, xfact=-1, yfact=1, origin=(0, 0))          # flipped over, as cad.py
    pins = [Point(RV.r_pins * cos(tau * k / RV.pins), RV.r_pins * sin(tau * k / RV.pins)) for k in range(RV.pins)]
    worst = (np.inf, -np.inf)
    for th_in in np.linspace(0, revs_in * tau, n):
        th_o, th_c = angles(th_in)
        for side, poly in ((1, disc), (-1, disc2)):
            c = side * disc_centre(th_o, th_c)
            g = affinity.translate(affinity.rotate(poly, th_o, origin=(0, 0), use_radians=True), c[0], c[1])
            gaps = [g.exterior.distance(pp) - RV.pin_r for pp in pins]
            inside = any(g.contains(pp) for pp in pins)
            if inside:
                return False, (-1.0, 0.0)
            worst = (min(worst[0], min(gaps)), max(worst[1], max(gaps)))
    ok = worst[0] > -0.005
    print(f"meshing check, {n} poses over {revs_in} input revolutions: disc-to-pin gap {worst[0] * 1e3:.1f} to "
          f"{worst[1] * 1e3:.1f} um (profile offset 30 um): {'ok' if ok else 'DISCS CUT THE PINS'}")
    return ok, worst


# ---- meshes: one per moving instance, in its build pose ------------------------------------------
STATIC = ["housing", "pin", "gear_cover", "motor_can", "stator", "end_windings", "input_bearing", "main_bearing"]


def prep_meshes():
    stamp = MESH / ".stamp"
    if stamp.exists() and stamp.stat().st_mtime > (OUT / "design.json").stat().st_mtime:
        return json.loads(stamp.read_text())
    from build123d import Compound, export_stl
    import cad
    MESH.mkdir(parents=True, exist_ok=True)
    print("exporting per-part meshes (phasing the planets first, a few minutes)...")
    cad.phase_by_trial()
    p = cad.assembly()
    names = {}

    def put(name, solids):
        export_stl(Compound(children=list(solids)), str(MESH / f"{name}.stl"), tolerance=0.05, angular_tolerance=0.3)
        names[name] = name

    for k in STATIC:
        put(k, p[k])
    for k in ("sun", "input_shaft", "rotor", "carrier_front", "carrier_rear", "journal_bearing"):
        put(k, p[k])
    for k in range(RV.n_crank):
        put(f"crank{k}", [p["crankshaft"][k]])
        put(f"planet{k}", [p["planet"][k]])
    put("disc1", [p["disc"][0]])
    put("disc2", [p["disc"][1]])
    n = RV.n_rollers
    i = 0
    for k in range(RV.n_crank):
        for side in (1, -1):
            put(f"rollers{k}{'a' if side > 0 else 'b'}", p["roller"][i:i + n])
            j = i // n
            put(f"eccring{k}{'a' if side > 0 else 'b'}", [p["ecc_ring"][j]])
            i += n
    stamp.write_text(json.dumps(names))
    return names


RGBA = dict(housing=".72 .74 .77 .18", pin=".88 .88 .9 1", gear_cover=".15 .16 .18 .18", motor_can=".15 .16 .18 .18",
            stator=".25 .27 .3 1", end_windings=".85 .45 .22 1", input_bearing=".85 .86 .88 1", main_bearing=".85 .86 .88 .25",
            sun=".8 .62 .35 1", input_shaft=".66 .68 .71 1", rotor=".2 .2 .22 1", carrier_front=".35 .45 .62 .22",
            carrier_rear=".35 .45 .62 .35", journal_bearing=".85 .86 .88 .5", crank=".55 .57 .6 1",
            planet=".8 .62 .35 1", eccring=".85 .86 .88 1", disc1=".82 .6 .28 1", disc2=".7 .45 .2 1", rollers=".93 .93 .95 1")


def colour(name):
    for k in sorted(RGBA, key=len, reverse=True):
        if name.startswith(k):
            return RGBA[k]
    return ".7 .7 .7 1"


MOVING = (["sun", "input_shaft", "rotor", "carrier_front", "carrier_rear", "journal_bearing", "disc1", "disc2"]
          + [f"crank{k}" for k in range(3)] + [f"planet{k}" for k in range(3)]
          + [f"rollers{k}{s}" for k in range(3) for s in "ab"] + [f"eccring{k}{s}" for k in range(3) for s in "ab"])


def make_xml(names):
    assets = "\n".join(f'<mesh name="{n}" file="{MESH / (n + ".stl")}" scale=".001 .001 .001"/>' for n in names)
    stat = "\n".join(f'<geom type="mesh" mesh="{n}" rgba="{colour(n)}" contype="0" conaffinity="0"/>' for n in STATIC)
    mov = "\n".join(f'<body name="{n}" mocap="true"><geom type="mesh" mesh="{n}" rgba="{colour(n)}" contype="0" conaffinity="0"/></body>'
                    for n in MOVING)
    # output marker: a tab on the carrier's output face; input marker: a flag on the motor end of the shaft
    z_end = 0.206
    return f"""<mujoco model="rv_actuator">
  <visual><global offwidth="1920" offheight="1080"/><quality shadowsize="4096"/>
    <headlight ambient=".35 .35 .35" diffuse=".6 .6 .6"/></visual>
  <asset>
    <texture name="sky" type="skybox" builtin="gradient" rgb1="1 1 1" rgb2=".80 .83 .87" width="512" height="512"/>
    {assets}
  </asset>
  <worldbody>
    <light pos=".3 -.4 .6" dir="-.4 .6 -.6" diffuse=".7 .7 .7"/>
    {stat}
    {mov}
    <body name="out_marker" mocap="true"><geom type="box" size=".012 .006 .003" pos=".062 0 {z_end + 0.003}" rgba="1 .25 .2 1" contype="0" conaffinity="0"/></body>
    <body name="in_marker" mocap="true"><geom type="box" size=".018 .003 .003" pos="0 .018 -.004" rgba=".2 .5 1 1" contype="0" conaffinity="0"/></body>
  </worldbody>
</mujoco>"""


def quat_z(a):
    return np.array([cos(a / 2), 0.0, 0.0, sin(a / 2)])


def apply(m, d, theta_in):
    th_o, th_c = angles(theta_in)
    poses = {"sun": pose_spin(theta_in), "input_shaft": pose_spin(theta_in), "rotor": pose_spin(theta_in),
             "in_marker": pose_spin(theta_in), "carrier_front": pose_spin(th_o), "carrier_rear": pose_spin(th_o),
             "journal_bearing": pose_spin(th_o), "out_marker": pose_spin(th_o),
             "disc1": pose_disc(1, th_o, th_c), "disc2": pose_disc(-1, th_o, th_c)}
    for k in range(3):
        poses[f"crank{k}"] = poses[f"planet{k}"] = pose_crank(k, th_o, th_c)
        poses[f"rollers{k}a"] = pose_rollers(k, 1, th_o, th_c)
        poses[f"rollers{k}b"] = pose_rollers(k, -1, th_o, th_c)
        poses[f"eccring{k}a"] = poses["disc1"]           # outer rings ride in the discs
        poses[f"eccring{k}b"] = poses["disc2"]
    for name, (a, t) in poses.items():
        mid = m.body_mocapid[m.body(name).id]
        d.mocap_pos[mid] = (t[0] * 1e-3, t[1] * 1e-3, 0.0)
        d.mocap_quat[mid] = quat_z(a)
    mujoco.mj_forward(m, d)
    return th_o, th_c


def build():
    names = prep_meshes()
    m = mujoco.MjModel.from_xml_string(make_xml(names))
    return m, mujoco.MjData(m)


_FONT = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf"


def _fonts():
    from PIL import ImageFont
    global FONT, FONT_B, FONT_S
    FONT_B = ImageFont.truetype(_FONT.replace("Sans.ttf", "Sans-Bold.ttf"), 26)
    FONT = ImageFont.truetype(_FONT, 22)
    FONT_S = ImageFont.truetype(_FONT, 17)


def video(seconds=12.0, fps=30, crank_rev_per_s=0.5):
    """crank_rev_per_s: crank turns relative to the carrier, per second of video. The output
    turns 1/40 of that, the input (i / 40) times it."""
    from PIL import Image, ImageDraw
    _fonts()
    ok, _ = check()
    if not ok:
        raise SystemExit("meshing check failed; not rendering")
    m, d = build()
    r = mujoco.Renderer(m, 720, 1280)
    cam = mujoco.MjvCamera()
    cam.lookat[:] = (0.0, 0.0, 0.165)
    cam.distance, cam.elevation = 0.40, -32
    w_in = crank_rev_per_s * tau * RV.z_planet / RV.z_sun * (1 - 1 / I) ** -1
    frames_dir = OUT / "sim_frames"
    frames_dir.mkdir(exist_ok=True)
    for f in frames_dir.glob("*.png"):
        f.unlink()
    n = int(seconds * fps)
    for i in range(n):
        t = i / fps
        th_in = w_in * t
        th_o, th_c = apply(m, d, th_in)
        cam.azimuth = 120 + 50 * sin(tau * t / seconds)
        r.update_scene(d, cam)
        im = Image.fromarray(r.render())
        dr = ImageDraw.Draw(im)
        dr.text((28, 22), f"RV actuator, {I:.1f}:1  (kinematic playback in MuJoCo)", fill=(20, 20, 20), font=FONT_B)
        dr.text((28, 62), f"input            {th_in / tau:7.2f} rev", fill=(30, 80, 200), font=FONT)
        dr.text((28, 90), f"crankshafts   {-th_c / tau:7.2f} rev  (relative to the carrier)", fill=(60, 60, 60), font=FONT)
        dr.text((28, 118), f"output (red tab) {np.degrees(th_o):6.2f} deg", fill=(200, 50, 40), font=FONT)
        dr.text((28, 680), "housing, carrier and covers transparent. Slowed so the crankshafts turn at "
                f"{crank_rev_per_s:g} rev/s (at rated speed: 9.75 rev/s)", fill=(70, 70, 70), font=FONT_S)
        im.save(frames_dir / f"f{i:04d}.png")
    mp4, gif = OUT / "rv_sim.mp4", OUT / "rv_sim.gif"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(fps), "-i", str(frames_dir / "f%04d.png"),
                    "-c:v", "libx264", "-pix_fmt", "yuv420p", "-crf", "20", str(mp4)], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-i", str(mp4), "-vf",
                    "fps=15,scale=720:-1:flags=lanczos,split[a][b];[a]palettegen=max_colors=128[p];[b][p]paletteuse",
                    str(gif)], check=True)
    print(f"wrote {mp4} ({mp4.stat().st_size / 1e6:.1f} MB) and {gif} ({gif.stat().st_size / 1e6:.1f} MB)")


def interactive(crank_rev_per_s=0.5):
    import time

    import mujoco.viewer
    m, d = build()
    w_in = crank_rev_per_s * tau * RV.z_planet / RV.z_sun * (1 - 1 / I) ** -1
    t0 = time.time()
    with mujoco.viewer.launch_passive(m, d) as v:
        while v.is_running():
            apply(m, d, w_in * (time.time() - t0))
            v.sync()
            time.sleep(1 / 60)


if __name__ == "__main__":
    {"interactive": interactive, "video": video, "check": lambda: check()}[_MODE]()
