"""
Blender side of explode.py: import the exploded GLB, light it, frame an orthographic camera
on it, render, and write where each named part group lands in the image (for the balloons).

    blender -b -P explode_blender.py -- in.glb out.png labels.json WIDTH HEIGHT SAMPLES
"""

import json
import math
import sys

import bpy
from bpy_extras.object_utils import world_to_camera_view
from mathutils import Vector

a = sys.argv[sys.argv.index("--") + 1:]
GLB, PNG, JS, W, H, SAMPLES = a[0], a[1], a[2], int(a[3]), int(a[4]), int(a[5])

bpy.ops.wm.read_factory_settings(use_empty=True)
sc = bpy.context.scene
bpy.ops.import_scene.gltf(filepath=GLB)
meshes = [o for o in sc.objects if o.type == "MESH"]
bpy.context.view_layer.update()


def corners(o):
    return [o.matrix_world @ Vector(c) for c in o.bound_box]


pts = [p for o in meshes for p in corners(o)]
lo = Vector((min(p.x for p in pts), min(p.y for p in pts), min(p.z for p in pts)))
hi = Vector((max(p.x for p in pts), max(p.y for p in pts), max(p.z for p in pts)))
center = (lo + hi) / 2
size = (hi - lo).length

# world: plain light grey, so the drawing reads like a document
world = bpy.data.worlds.new("w")
sc.world = world
world.use_nodes = True
bg = world.node_tree.nodes["Background"]
bg.inputs[0].default_value = (0.96, 0.96, 0.97, 1)
bg.inputs[1].default_value = 0.9


def area(name, loc, energy, s):
    L = bpy.data.lights.new(name, "AREA")
    L.energy, L.size = energy, s
    ob = bpy.data.objects.new(name, L)
    ob.location = loc
    sc.collection.objects.link(ob)
    d = center - ob.location
    ob.rotation_euler = d.to_track_quat("-Z", "Y").to_euler()


area("key", center + Vector((size * 0.3, -size * 0.6, size * 0.7)), size * size * 2.5, size * 0.6)
area("fill", center + Vector((-size * 0.5, -size * 0.4, size * 0.2)), size * size * 1.0, size * 0.8)

# camera: orthographic, looking from front-left-above so the disc faces show (axis along x)
cam = bpy.data.objects.new("cam", bpy.data.cameras.new("cam"))
sc.collection.objects.link(cam)
sc.camera = cam
cam.data.type = "ORTHO"
cam.data.clip_start, cam.data.clip_end = size * 0.01, size * 10   # mm read as m: the model is hundreds of units long
d = Vector((-0.30, -0.86, 0.40)).normalized()
cam.location = center + d * size * 2.0
cam.rotation_euler = (center - cam.location).to_track_quat("-Z", "Y").to_euler()
sc.render.resolution_x, sc.render.resolution_y = W, H
bpy.context.view_layer.update()
# fit the ortho scale to the projected extent
cam.data.ortho_scale = size
bpy.context.view_layer.update()
uv = [world_to_camera_view(sc, cam, p) for p in pts]
span_x = max(u.x for u in uv) - min(u.x for u in uv)
span_y = (max(u.y for u in uv) - min(u.y for u in uv)) * H / W
cam.data.ortho_scale = size * max(span_x, span_y) * 1.08
bpy.context.view_layer.update()
uv = [world_to_camera_view(sc, cam, p) for p in pts]
mid_u = (max(u.x for u in uv) + min(u.x for u in uv)) / 2
mid_v = (max(u.y for u in uv) + min(u.y for u in uv)) / 2
cam.data.shift_x, cam.data.shift_y = mid_u - 0.5, (mid_v - 0.5) * H / W

sc.render.engine = "CYCLES"
sc.cycles.samples = SAMPLES
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
        break
sc.view_settings.view_transform = "Standard"
sc.render.filepath = PNG
bpy.ops.render.render(write_still=True)

# label anchors: each group's projected centre and its top-most projected point
bpy.context.view_layer.update()
out = {}
for o in meshes:
    name = o.name.split(".")[0]
    cs = corners(o)
    c = sum(cs, Vector()) / 8
    uc = world_to_camera_view(sc, cam, c)
    # anchors on the part itself: of the mesh vertices facing the camera side, the one whose
    # projection is nearest the projected centre, and the highest and lowest projected ones
    mw = o.matrix_world
    vs = [world_to_camera_view(sc, cam, mw @ v.co) for v in o.data.vertices[::max(1, len(o.data.vertices) // 4000)]]
    vs = [u for u in vs if u.z > 0]
    near = min(vs, key=lambda u: (u.x - uc.x) ** 2 + ((u.y - uc.y) * H / W) ** 2)
    top = max(vs, key=lambda u: u.y)
    bot = min(vs, key=lambda u: u.y)
    px = lambda u: (u.x * W, (1 - u.y) * H)
    out[name] = dict(center=px(near), top=px(top), bottom=px(bot))
json.dump(out, open(JS, "w"), indent=1)
print("WROTE", PNG, JS)
