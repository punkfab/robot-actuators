"""
Blender side of assemble_anim.py. Keys every part's pose per frame (exploded -> assembled,
then the reducer running), fades the outer parts, moves an orthographic camera, renders frames.

    blender -b -P assemble_blender.py -- parts.glb meta.json frames_dir
"""

import json
import math
import sys

import bpy
from mathutils import Matrix, Vector

GLB, JS, FRAMES = sys.argv[sys.argv.index("--") + 1:][:3]
M = json.load(open(JS))
W, H = M["res"]

bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
bpy.ops.import_scene.gltf(filepath=GLB)
bpy.context.view_layer.update()
objs = {o.name.split(".")[0]: o for o in sc.objects if o.type == "MESH"}
info = {d["name"]: d for d in M["objects"]}

# the importer maps glTF (x, y, z) to Blender (x, -z, y); C undoes that, L lays the actuator
# along Blender's X with CAD -x up
C = Matrix(((1, 0, 0, 0), (0, 0, -1, 0), (0, 1, 0, 0), (0, 0, 0, 1)))
L = Matrix(((0, 0, 1, 0), (0, 1, 0, 0), (-1, 0, 0, 0), (0, 0, 0, 1)))
M0 = {n: o.matrix_world.copy() for n, o in objs.items()}
for o in objs.values():
    o.parent = None
    o.rotation_mode = "QUATERNION"


def T(angle, tx, ty, tz):
    """CAD-space rigid transform: rotate about +z by angle, then translate"""
    return Matrix.Translation((tx, ty, tz)) @ Matrix.Rotation(angle, 4, "Z")


def about(point, angle):
    """rotation about +z through a point (x, y)"""
    px, py = point
    return Matrix.Translation((px, py, 0)) @ Matrix.Rotation(angle, 4, "Z") @ Matrix.Translation((-px, -py, 0))


def crank_xy(k):
    a = math.radians(M["crank0"]) + 2 * math.pi * k / M["n_crank"]
    return M["r_crank"] * math.cos(a), M["r_crank"] * math.sin(a)


def kinematic(role, k, th_in):
    """CAD-space pose of a part at input angle th_in (the reducer's exact kinematics)"""
    i, e = M["ratio"], M["e"]
    th_o = th_in / i
    th_c = -(M["z_sun"] / M["z_planet"]) * (th_in - th_o)
    R_o = Matrix.Rotation(th_o, 4, "Z")
    if role == "fixed":
        return Matrix.Identity(4)
    if role == "input":
        return Matrix.Rotation(th_in, 4, "Z")
    if role == "carrier":
        return R_o
    if role == "crank":
        return R_o @ about(crank_xy(k), th_c)
    if role in ("disc1", "disc2"):
        side = 1 if role == "disc1" else -1
        phi = th_o + th_c
        c = (side * e * math.cos(phi), side * e * math.sin(phi))
        # rotate about the built centre (side*e, 0), then move that centre to c
        return Matrix.Translation((c[0] - side * e, c[1], 0)) @ about((side * e, 0.0), th_o)
    raise ValueError(role)


def ease(x):
    x = min(1.0, max(0.0, x))
    return x * x * (3 - 2 * x)


t0, stg, mv, act2, run = M["t0"], M["stagger"], M["move"], M["act2"], M["run"]
n_frames = act2 + run
w_in = M["crank_revs"] * 2 * math.pi / run / (M["z_sun"] / M["z_planet"]) / (1 - 1 / M["ratio"])


def pose(name, f):
    d = info[name]
    g = d["group"]
    start = t0 + g * stg
    s = 1.0 - ease((f - start) / mv)                       # 1 exploded .. 0 home
    dz = M["dz"][str(g)] * s
    th_in = max(0, f - act2) * w_in
    return L @ Matrix.Translation((0, 0, dz)) @ kinematic(d["role"], d["k"], th_in) @ C.inverted() @ M0[name]


# fading materials: housing, cover, main bearings go to glass during act 2
for name, o in objs.items():
    if not info[name]["fade"]:
        continue
    mat = o.active_material.copy()
    o.active_material = mat
    bsdf = next(n for n in mat.node_tree.nodes if n.type == "BSDF_PRINCIPLED")
    a = bsdf.inputs["Alpha"]
    for f, v in ((act2 - 5, 1.0), (act2 + 25, 0.18)):
        a.default_value = v
        a.keyframe_insert("default_value", frame=f)

for f in range(n_frames):
    for name, o in objs.items():
        o.matrix_world = pose(name, f)
        o.keyframe_insert("location", frame=f)
        o.keyframe_insert("rotation_quaternion", frame=f)


# camera: orthographic, fit to the exploded layout at the start, to the assembly at the end of
# act 1, then a slow orbit about the actuator axis (Blender X) during act 2
def world_points(f):
    pts = []
    for name, o in objs.items():
        m = pose(name, f)
        pts += [m @ Vector(c) for c in o.bound_box]
    return pts


cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
sc.collection.objects.link(cam)
sc.camera = cam
cam.data.type = "ORTHO"


def fit(pts, direction):
    lo = Vector((min(p.x for p in pts), min(p.y for p in pts), min(p.z for p in pts)))
    hi = Vector((max(p.x for p in pts), max(p.y for p in pts), max(p.z for p in pts)))
    centre = (lo + hi) / 2
    size = (hi - lo).length
    q = (-direction).to_track_quat("-Z", "Y")
    right, up = q @ Vector((1, 0, 0)), q @ Vector((0, 1, 0))
    xs = [(p - centre).dot(right) for p in pts]
    ys = [(p - centre).dot(up) for p in pts]
    cx, cy = (max(xs) + min(xs)) / 2, (max(ys) + min(ys)) / 2
    scale = max(max(xs) - min(xs), (max(ys) - min(ys)) * W / H) * 1.10
    centre = centre + right * cx + up * cy
    return centre, scale, size, q


IN_SIDE = Vector((-0.34, -0.84, 0.42)).normalized()     # act 1: from the motor end, as explode.py
OUT_SIDE = Vector((0.62, -0.62, 0.48)).normalized()     # act 2: from the output end, motor behind


def view_dir(orbit, swing=0.0):
    """unit vector from the target toward the camera"""
    base = IN_SIDE.lerp(OUT_SIDE, swing).normalized()
    return Matrix.Rotation(orbit, 3, "X") @ base


c_ex, s_ex, size, q = fit(world_points(0), view_dir(0))
c_as, s_as, _, _ = fit(world_points(act2), view_dir(0, 1.0))
cam.data.clip_start, cam.data.clip_end = size * 0.01, size * 20
# act 1: fit the camera to the parts as they are on each frame (exploded and arriving), smoothed;
# act 2: hold the assembled fit and orbit slowly about the actuator axis
fits = [fit(world_points(f), view_dir(0))[:2] for f in range(0, act2, 3)]
centre_s, scale_s = fits[0]
for f in range(n_frames):
    if f < act2:
        c_f, s_f = fits[min(f // 3, len(fits) - 1)]
        a = 0.12
        centre_s = centre_s.lerp(c_f, a)
        scale_s = scale_s + (s_f - scale_s) * a
        centre, scale, orbit, swing = centre_s, scale_s, 0.0, 0.0
        if f > act2 - 40:                                   # settle onto the act-2 framing, swinging round
            u = ease((f - (act2 - 40)) / 40)
            centre, scale, swing = centre.lerp(c_as, u), scale + (s_as * 1.15 - scale) * u, u
    else:
        u = (f - act2) / run
        orbit = math.radians(-50) * ease(u)
        centre, scale, swing = c_as, s_as * 1.15, 1.0
    d = view_dir(orbit, swing)                                     # from the target toward the camera
    cam.location = centre + d * size * 2.5
    cam.rotation_mode = "QUATERNION"
    cam.rotation_quaternion = (-d).to_track_quat("-Z", "Y")
    cam.data.ortho_scale = scale
    cam.keyframe_insert("location", frame=f)
    cam.keyframe_insert("rotation_quaternion", frame=f)
    cam.data.keyframe_insert("ortho_scale", frame=f)

# lights follow the camera's orbit loosely: a big soft key and a fill
world = bpy.data.worlds.new("w")
sc.world = world
world.use_nodes = True
bg = world.node_tree.nodes["Background"]
bg.inputs[0].default_value = (0.95, 0.95, 0.96, 1)
bg.inputs[1].default_value = 0.85
for name, off, e, sz in (("key", Vector((0.3, -0.7, 0.8)), 3.0, 0.7), ("fill", Vector((-0.6, -0.3, 0.2)), 1.2, 0.9),
                         ("rim", Vector((0.2, 0.8, 0.5)), 1.5, 0.5)):
    Lt = bpy.data.lights.new(name, "AREA")
    Lt.energy, Lt.size = size * size * e, size * sz
    ob = bpy.data.objects.new(name, Lt)
    ob.location = c_as + off * size * 1.5
    ob.rotation_euler = (c_as - ob.location).to_track_quat("-Z", "Y").to_euler()
    sc.collection.objects.link(ob)

sc.render.engine = "CYCLES"
sc.cycles.samples = M["samples"]
sc.cycles.use_denoising = True
prefs = bpy.context.preferences.addons["cycles"].preferences
for backend in ("OPTIX", "CUDA", "HIP", "ONEAPI", "METAL"):
    prefs.compute_device_type = backend
    try:
        prefs.refresh_devices()
    except Exception:
        continue
    if any(dv.type == backend for dv in prefs.devices):
        for dv in prefs.devices:
            dv.use = dv.type == backend
        sc.cycles.device = "GPU"
        print("CYCLES device:", backend)
        break
sc.view_settings.view_transform = "Standard"
sc.render.resolution_x, sc.render.resolution_y = W, H
sc.frame_start, sc.frame_end = 0, n_frames - 1
sc.render.fps = M["fps"]
sc.render.filepath = FRAMES + "/f"
sc.render.image_settings.file_format = "PNG"
bpy.ops.render.render(animation=True)
print("RENDERED", n_frames)
