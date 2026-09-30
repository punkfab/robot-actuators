"""
Re-run a software-mfg gate with the STS3215 modelled by BAM instead of the ideal servo.

software-mfg/scripts/workcell_check.py solves IK to a point 60 mm above the work datum,
drives the SO-101 there for 4 s under the MJCF's stock position actuators, and passes
if the jaw arrives within 20 mm. This runs the same gate unchanged except for the
actuator model (and optional payloads), without touching software-mfg:

    python sts3215/crosscheck.py

It's the template for re-checking any punkfab gate that steps an SO-101: build the
spec the gate builds, bamify() the arm actuators, drive the same targets with StsServos.
"""

import sys
from pathlib import Path

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from arm import JOINTS, StsServos, bamify  # noqa: E402

SWMFG = Path.home() / "sandbox/punkfab/software-mfg"
sys.path.insert(0, str(SWMFG / "sim"))
import ik  # noqa: E402
from workcell import DATUM_POS, build_model, build_spec  # noqa: E402

TOL = 0.02          # workcell_check.py's pass threshold (m)
STEPS = 2000        # its drive duration (2000 × 2 ms)


def reach(kind, payload_kg=0.0, **servo_kw):
    spec = build_spec()
    if payload_kg:
        jaw = next(b for b in spec.bodies if b.name == "moving_jaw_so101_v1")
        g = jaw.parent.add_body(name="payload", pos=[0, 0, -0.09])
        g.add_geom(type=mujoco.mjtGeom.mjGEOM_SPHERE, size=[0.015, 0, 0], mass=payload_kg,
                   contype=0, conaffinity=0)
    if kind != "stock":
        bamify(spec, JOINTS)
    m = spec.compile()
    tip = m.body("moving_jaw_so101_v1").id
    target = DATUM_POS + np.array([0.0, 0.0, 0.06])
    q_arm, _ = ik.solve_ik(build_model(), mujoco.MjData(build_model()), tip, target)
    d = mujoco.MjData(m)
    mujoco.mj_forward(m, d)
    if kind == "stock":
        d.ctrl[:ik.ARM_DOFS] = q_arm
        mujoco.mj_step(m, d, STEPS)
    else:
        s = StsServos(m, d, JOINTS, kind, **servo_kw)
        for j, q in zip(JOINTS, list(q_arm) + [d.qpos[m.jnt_qposadr[m.joint("gripper").id]]]):
            s.set(j, q)
        for _ in range(STEPS):
            s.update()
            mujoco.mj_step(m, d)
    return float(np.linalg.norm(d.xpos[tip] - target))


def main():
    rng = np.random.default_rng(1)
    print("software-mfg workcell gate: jaw to 60 mm above the datum, pass < 20 mm\n")
    print(f"  {'payload':>7s}  {'stock':>7s}  {'BAM M6':>7s}  {'M6 band (friction ×0.7–1.3, 6.8–8 V)':>38s}")
    for pay in (0.0, 0.1, 0.2, 0.3):
        e_s, e_b = reach("stock", pay), reach("m6", pay)
        band = [reach("m6", pay, friction_scale=rng.uniform(0.7, 1.3), vin=rng.uniform(6.8, 8.0))
                for _ in range(12)]
        fails = sum(b > TOL for b in band)
        print(f"  {pay*1000:5.0f} g  {e_s*1000:6.1f}mm {e_b*1000:6.1f}mm   "
              f"{min(band)*1000:5.1f}–{max(band)*1000:5.1f} mm, {fails}/12 would FAIL"
              f"{'   <- stock PASS, M6 FAIL' if e_s < TOL <= e_b else ''}")


if __name__ == "__main__":
    main()
