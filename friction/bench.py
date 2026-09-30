"""
Pendulum test bench in MuJoCo: checks that FrictionUpdater makes the sim obey the
friction model, for several models (M1, our M5 prior, BAM-fitted STS3215 M4–M6).

  1. Static drive/backdrive window. At each gravity load τe, bisect in the SIM for
     the motor torques where the arm starts to lift (drive) and to fall (backdrive),
     and compare with Friction.hold_window() (the model's prediction).
  2. Moving friction. With a brake and viscous damping, the arm settles at a terminal
     velocity that the model predicts in closed form. Covers drive, backdrive, Stribeck.

  python friction/bench.py          # report + out/drive_backdrive.png
"""

import sys
from pathlib import Path

import mujoco
import numpy as np
from scipy.optimize import brentq

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from models import Friction  # noqa: E402
from mujoco_friction import FrictionUpdater  # noqa: E402

OUT = HERE / "out"
L = 0.10         # m arm length
DT = 0.0005      # s

XML = f"""
<mujoco model="friction_bench">
  <option timestep="{DT}" gravity="0 0 -9.81" integrator="implicitfast"/>
  <worldbody>
    <body name="arm">
      <joint name="j" type="hinge" axis="0 -1 0" armature="0.002"/>
      <geom type="capsule" fromto="0 0 0 {L} 0 0" size="0.004" mass="0.001"/>
      <body name="mass" pos="{L} 0 0">
        <geom type="sphere" size="0.01" mass="0.5"/>
      </body>
    </body>
  </worldbody>
  <actuator><motor name="m" joint="j" gear="1"/></actuator>
</mujoco>
"""

# Each model is exercised at its own sensible torque scale (a hobby servo is not a 40:1 drive).
MODELS = {
    "M1 Coulomb-viscous": (Friction.m1(kc=0.08, kv=0.01), 0.5),
    "M5 prior (η 83%)":   (Friction.from_efficiency(0.83, drag=0.08, kv=0.01), 0.5),
    "STS3215 M4 (BAM)":   (Friction.from_bam("sts3215", "m4")[0], 0.5),
    "STS3215 M5 (BAM)":   (Friction.from_bam("sts3215", "m5")[0], 0.5),
    "STS3215 M6 (BAM)":   (Friction.from_bam("sts3215", "m6")[0], 0.5),
}


def make(fr: Friction, mass: float):
    m = mujoco.MjModel.from_xml_string(XML)
    m.body_mass[m.body("mass").id] = mass
    d = mujoco.MjData(m)
    up = FrictionUpdater(m, "j", fr)
    return m, d, up


def moves(fr, mass, tau_m, t=0.3):
    """+1 lifts, −1 falls, 0 holds — from rest at horizontal with torque tau_m."""
    m, d, up = make(fr, mass)
    d.ctrl[0] = tau_m
    up.prime(d)
    up.step(d, int(t / DT))
    q = d.qpos[0]
    return 1 if q > 2e-3 else (-1 if q < -2e-3 else 0)


def sim_window(fr, mass):
    """Bisect the sim for the backdrive / drive thresholds at horizontal."""
    e = mass * 9.81 * L
    lo_a, lo_b = 0.0, e            # falls at lo_a, holds (or lifts) at lo_b
    if moves(fr, mass, 0.0) >= 0:  # holds unpowered: backdrive threshold is 0
        lo = 0.0
    else:
        for _ in range(14):
            mid = 0.5 * (lo_a + lo_b)
            lo_a, lo_b = (mid, lo_b) if moves(fr, mass, mid) < 0 else (lo_a, mid)
        lo = 0.5 * (lo_a + lo_b)
    hi_a, hi_b = e, 4 * e + 1.0    # holds at hi_a, lifts at hi_b
    for _ in range(14):
        mid = 0.5 * (hi_a + hi_b)
        hi_a, hi_b = (hi_a, mid) if moves(fr, mass, mid) > 0 else (mid, hi_b)
    return lo, 0.5 * (hi_a + hi_b)


def terminal_velocity(fr, tau_m, brake):
    """Gravity off, brake torque −brake, motor tau_m: sim vs closed-form terminal speed."""
    m, d, up = make(fr, 0.5)
    m.opt.gravity[:] = 0
    d.ctrl[0] = tau_m
    d.qfrc_applied[0] = -brake
    t = 10 * (0.002 + 0.5 * L**2) / fr.kv      # 10 time constants J/kv
    up.prime(d)
    up.step(d, int(t / DT))
    w_sim = d.qvel[0]

    def resid(w):   # kv·w = τm − brake − sign(w)·budget(τm, −brake, w)
        return fr.kv * w - (tau_m - brake - np.sign(w) * fr.budget(tau_m, -brake, w))
    try:
        w_model = brentq(resid, 1e-9, 500) if resid(1e-9) < 0 else brentq(resid, -500, -1e-9)
    except ValueError:
        w_model = 0.0
    return w_sim, w_model


def main():
    OUT.mkdir(exist_ok=True)
    print("=== 1. static drive/backdrive window at horizontal (sim bisection vs model) ===")
    print(f"  {'model':20s} {'τe':>6s}   {'backdrive below':>16s}   {'drive above':>14s}")
    diagram = {}
    worst = 0.0
    for name, (fr, _) in MODELS.items():
        pts = []
        for mass in (0.25, 0.5, 1.0):
            e = mass * 9.81 * L
            lo_s, hi_s = sim_window(fr, mass)
            lo_m, hi_m = fr.hold_window(e)
            pts.append((e, lo_s, hi_s, lo_m, hi_m))
            err = max(abs(lo_s - lo_m), abs(hi_s - hi_m)) / e
            worst = max(worst, err)
            print(f"  {name:20s} {e:6.3f}   sim {lo_s:5.3f} mdl {lo_m:5.3f}   "
                  f"sim {hi_s:5.3f} mdl {hi_m:5.3f}")
        diagram[name] = (fr, pts)
    print(f"  worst mismatch {worst:.1%} of the load  "
          f"{'✓' if worst < 0.03 else '✗ (check solver softness)'}\n")

    print("=== 2. moving friction: terminal velocity (gravity off, brake 0.2 N·m) ===")
    print(f"  {'model':20s} {'τm':>5s} {'state':>9s} {'ω sim':>8s} {'ω model':>8s}")
    worst_w = 0.0
    for name, (fr, _) in MODELS.items():
        if fr.kv <= 0:
            continue
        for tau_m, state in ((0.6, "drive"), (0.02, "backdrive")):
            w_s, w_m = terminal_velocity(fr, tau_m, 0.2)
            worst_w = max(worst_w, abs(w_s - w_m) / max(abs(w_m), 1.0))
            print(f"  {name:20s} {tau_m:5.2f} {state:>9s} {w_s:8.2f} {w_m:8.2f}")
    print(f"  worst mismatch {worst_w:.1%}  {'✓' if worst_w < 0.03 else '✗'}\n")

    plot(diagram)
    print(f"wrote {OUT/'drive_backdrive.png'}")


def plot(diagram):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    show = ["M1 Coulomb-viscous", "M5 prior (η 83%)", "STS3215 M6 (BAM)"]
    fig, axes = plt.subplots(1, len(show), figsize=(4.2 * len(show), 4), sharey=True)
    for ax, name in zip(axes, show):
        fr, pts = diagram[name]
        e = np.linspace(0, 1.0, 60)
        win = np.array([fr.hold_window(x) for x in e])
        ax.fill_betweenx(e, win[:, 0], win[:, 1], color="0.85", label="holds (static)")
        ax.plot(win[:, 1], e, "C0", label="drive threshold (model)")
        ax.plot(win[:, 0], e, "C3", label="backdrive threshold (model)")
        ax.plot(e, e, "k--", lw=0.8, label="τm = τe (no friction)")
        P = np.array(pts)
        ax.plot(P[:, 2], P[:, 0], "o", mfc="none", color="C0", label="sim")
        ax.plot(P[:, 1], P[:, 0], "o", mfc="none", color="C3")
        ax.set_title(f"{name}\n{fr.describe().split('|')[1].strip()}", fontsize=9)
        ax.set_xlabel("motor torque τm at output (N·m)")
        ax.set_xlim(0, 1.6)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("load τe (N·m)")
    axes[0].legend(fontsize=7, loc="lower right")
    fig.suptitle("Drive/backdrive diagrams: where the load holds, lifts, or falls", fontsize=11)
    fig.tight_layout()
    fig.savefig(OUT / "drive_backdrive.png", dpi=130)


if __name__ == "__main__":
    main()
