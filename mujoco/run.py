"""
Drive the cycloidal actuator in MuJoCo and print a sizing report.

  python mujoco/run.py            # headless physics + report
  python mujoco/run.py --m1       # legacy friction (gear = η∞, Coulomb drag) for A/B
  python mujoco/run.py --view     # interactive viewer (needs a display)

Injects armature / gear / friction / payload from actuator.py so the sim always
matches the current CAD Params and motor constants.
"""

import sys
from math import pi
from pathlib import Path

import mujoco
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "cycloidal"))
sys.path.insert(0, str(HERE.parent / "friction"))
from actuator import MotorSpec, ActuatorSpec, Params  # noqa: E402
from efficiency import predict_params  # noqa: E402  (Layer B -> Layer A)
from models import Friction  # noqa: E402
from mujoco_friction import FrictionUpdater  # noqa: E402

# --- losses + test payload (knobs) ------------------------------------------- #
# Gearbox friction is an extended (BAM-style) model, updated every step by
# friction/mujoco_friction.py. The gear carries the IDEAL Kt·N and the friction
# budget carries every loss:
#     budget = Kc + |Km·τm − Ke·τe|      (M5 prior from the Layer-B η∞ + drag)
# Kc is chosen so the drive curve η(T)=η∞·T/(T+drag) is unchanged, but holding
# and backdriving now see load-dependent friction too (see friction/README.md).
# ETA_BACK = None -> symmetric (η_back = η∞), the unbiased prior until measured.
#
# --m1 restores the old model: gear = η∞·Kt·N, frictionloss = drag. That is
# right when driving, but makes the motor pay τe/η just to HOLD, and backdrives
# with only the no-load drag resisting.
LEGACY_M1       = "--m1" in sys.argv
_SPEC0          = ActuatorSpec.from_motor(MotorSpec(), Params())
FRICTIONLOSS    = _SPEC0.drag_out   # N·m no-load drag at the output
DAMPING         = 0.0015   # N·m·s/rad viscous
ETA_BACK        = None     # backdrive efficiency; None = symmetric prior
PAYLOAD_KG      = 0.10     # mass at the end of the 150 mm test arm
ARM_LEN         = 0.15     # m  (must match testbench.xml)


def gearbox_friction(spec, legacy=LEGACY_M1) -> Friction:
    if legacy:
        return Friction.m1(kc=spec.drag_out, kv=DAMPING)
    return Friction.from_efficiency(spec.eta_inf, spec.drag_out, kv=DAMPING, eta_back=ETA_BACK)


def load(legacy=LEGACY_M1):
    motor = MotorSpec()
    spec = ActuatorSpec.from_motor(motor, Params())
    model = mujoco.MjModel.from_xml_path(str(HERE / "testbench.xml"))

    # inject derived physics onto the actuated joint + actuator
    jid = model.joint("joint_out").id
    dof = model.jnt_dofadr[jid]
    model.dof_armature[dof] = spec.reflected_inertia

    aid = model.actuator("drive").id
    # legacy: gear carries η∞; extended: gear is ideal and the friction carries the loss
    model.actuator_gear[aid, 0] = spec.torque_per_amp if legacy else motor.kt * spec.ratio
    model.actuator_ctrlrange[aid] = [-MotorSpec().max_current, MotorSpec().max_current]
    fr = FrictionUpdater(model, "joint_out", gearbox_friction(spec, legacy))

    # set the test payload
    pid = model.body("payload").id
    model.body_mass[pid] = PAYLOAD_KG

    return model, spec, motor, aid, fr


def new_data(model, fr):
    data = mujoco.MjData(model)
    fr.prime(data)
    return data


def settle_angle(model, data, ctrl, t, aid, fr):
    data.ctrl[aid] = ctrl
    fr.step(data, int(t / model.opt.timestep))
    return data.qpos[0], data.qvel[0]


def drive_voltage(model, data, throttle, t, aid, motor, N, dof, fr):
    """Step with a throttle (voltage) command; current is back-EMF limited each step."""
    n = int(t / model.opt.timestep)
    for _ in range(n):
        data.ctrl[aid] = motor.current_at(data.qvel[dof] * N, throttle)
        fr.step(data)
    return data.qpos[dof], data.qvel[dof]


def hold_window_sim(model, aid, fr, t=0.3):
    """Bisect the sim for the currents that just hold the payload at horizontal:
    below lo it backdrives (falls), above hi it drives (lifts)."""
    d0 = new_data(model, fr)
    up = -np.sign(-d0.qfrc_bias[0])           # "lift" = against gravity

    def moves(I):                            # +1 lifts, −1 falls, 0 holds
        d = new_data(model, fr)
        q, _ = settle_angle(model, d, up * I, t, aid, fr)
        q *= up
        return 1 if q > 2e-3 else (-1 if q < -2e-3 else 0)
    Imax = MotorSpec().max_current
    lo = 0.0
    if moves(0.0) < 0:
        a, b = 0.0, Imax
        for _ in range(16):
            mid = 0.5 * (a + b)
            a, b = (mid, b) if moves(mid) < 0 else (a, mid)
        lo = 0.5 * (a + b)
    a, b = lo, Imax
    for _ in range(16):
        mid = 0.5 * (a + b)
        a, b = (a, mid) if moves(mid) > 0 else (mid, b)
    return lo, 0.5 * (a + b)


def main():
    model, spec, motor, aid, fr = load()
    spec.report(motor)
    print(f"gearbox friction: {fr.f.describe()}"
          f"{'   [--m1 legacy: gear carries η∞]' if LEGACY_M1 else ''}\n")
    dof = model.jnt_dofadr[model.joint("joint_out").id]
    Kt, N = motor.kt, spec.ratio
    gear = model.actuator_gear[aid, 0]

    # --- Scenario A: free angular acceleration (gravity off) -> verify inertia --
    model.opt.gravity[:] = 0
    data = new_data(model, fr)
    _, v = settle_angle(model, data, MotorSpec().max_current, 0.05, aid, fr)
    # net torque = actuator torque minus the friction budget it slides against
    tau_act = gear * MotorSpec().max_current
    net_t = tau_act - float(fr.f.total(tau_act, 0.0, v / 2))   # mean speed over the ramp
    I_eff = net_t / (v / 0.05) if v > 1e-6 else float("nan")
    print("--- Scenario A: free accel (gravity off, peak current) ---")
    print(f"  output reached {v:6.1f} rad/s in 50 ms  ->  effective inertia {I_eff*1e4:.2f}e-4 kg·m²")
    print(f"  (reflected motor inertia alone = {spec.reflected_inertia*1e4:.2f}e-4; rest is arm+payload)\n")

    # --- Scenario B: lift the test arm from horizontal (gravity on, full throttle) --
    model.opt.gravity[:] = (0, 0, -9.81)
    data = new_data(model, fr)
    grav_torque = data.qfrc_bias[0]   # torque needed to hold at horizontal
    q, w = drive_voltage(model, data, 1.0, 0.6, aid, motor, N, dof, fr)
    lifted = np.degrees(q)
    print("--- Scenario B: lift 100 g @ 150 mm from horizontal (full throttle) ---")
    print(f"  gravity hold torque needed = {abs(grav_torque):.3f} N·m   (peak avail {spec.peak_torque:.3f})")
    print(f"  arm swung to {lifted:+.0f}° at {abs(w)*60/2/pi:.0f} rpm  ->  "
          f"{'LIFTS ✓' if lifted > 60 else 'STALLS ✗'}  (back-EMF now caps the speed)\n")

    # --- Scenario C: backdrivability (no power, gravity on) --------------------
    data = new_data(model, fr)
    q, _ = settle_angle(model, data, 0.0, 0.6, aid, fr)
    drop = np.degrees(q)
    # unpowered: friction budget at τm = 0 is Kc + Ke·|τe| (load-dependent)
    f_hold = float(fr.f.budget(0.0, grav_torque, 0.0))
    held = abs(grav_torque) <= f_hold
    print("--- Scenario C: unpowered hold (backdrive) ---")
    print(f"  friction budget {f_hold:.3f} N·m vs gravity {abs(grav_torque):.3f} N·m  ->  "
          f"{'self-holds' if held else f'backdrives, falls to {drop:+.0f}°'}")
    if LEGACY_M1:
        print("  (legacy: only the no-load drag resists backdriving, whatever the load)\n")
    else:
        print(f"  (unpowered it self-holds up to {fr.f.kc/max(1e-9, 1-fr.f.ke):.3f} N·m; "
              f"η_back {fr.f.eta_back:.0%} is a prior until measured)\n")

    # --- Scenario D: efficiency vs load, measured from sim power balance --------
    # For each output torque T, command the current that delivers it, apply T as a
    # brake, run to steady speed, and read η = P_out/P_in from the sim. Demonstrates
    # the torque-based model (η rising with load) is physically in the sim.
    print("--- Scenario D: efficiency vs load (measured in sim) ---")
    print(f"  {'T_out':>7s} {'I (A)':>6s} {'out rpm':>8s} {'η meas':>7s} {'η model':>8s}")
    model.opt.gravity[:] = 0
    for T in (0.05, 0.10, 0.20, 0.40, spec.peak_torque):
        # current whose torque just drives T through the friction, +1% so it creeps forward
        I = fr.f.hold_window(T, dq=0.3)[1] / gear * 1.01
        if I > motor.max_current + 1e-9:
            print(f"  {T:7.3f}  {'>13':>5s}   {'—':>7s}    —      {spec.eta_at(T)*100:5.0f}%  (exceeds 13 A)")
            continue
        data = new_data(model, fr)
        data.qfrc_applied[dof] = -T              # brake opposing the (positive) motion
        settle_angle(model, data, I, 0.8, aid, fr)
        w = data.qvel[dof]
        p_out = T * w
        p_in = Kt * I * (N * w)
        eta_meas = p_out / p_in if p_in > 1e-9 else 0.0
        print(f"  {T:7.3f} {I:6.1f} {w*60/(2*pi):8.0f} {eta_meas*100:6.0f}% {spec.eta_at(T)*100:7.0f}%")
    print("  (η meas tracks η model -> the load-dependent curve emerges from the friction budget)\n")

    # --- Scenario E: torque-speed curve traced from a real spin-up (back-EMF) ----
    # Full throttle, no load. As the output accelerates, back-EMF cuts the current,
    # so the delivered torque is flat then droops to zero near no-load speed.
    print("--- Scenario E: torque-speed curve (back-EMF, traced in sim spin-up) ---")
    model.opt.gravity[:] = 0
    pid = model.body("payload").id
    model.body_mass[pid] = 0.001                  # light, so it reaches high speed quickly
    data = new_data(model, fr)
    # sample the curve at fractions of the (ratio-dependent) no-load speed
    nl_rpm = spec.no_load_speed * 60 / (2 * pi)
    targets, got = [round(f * nl_rpm) for f in (0.0, 0.25, 0.5, 0.7, 0.85, 0.95)], {}
    for _ in range(int(4.0 / model.opt.timestep)):
        rpm = data.qvel[dof] * 60 / (2 * pi)
        I = motor.current_at(data.qvel[dof] * N, 1.0)
        tau = fr.f.load_torque(gear * I, data.qvel[dof])   # load it could drive here
        for tgt in targets:
            if tgt not in got and rpm >= tgt:
                got[tgt] = tau
        data.ctrl[aid] = I
        fr.step(data)
    model.body_mass[pid] = PAYLOAD_KG             # restore
    print(f"  {'rpm':>6s} {'τ meas':>8s} {'τ model':>8s}")
    for tgt in targets:
        if tgt in got:
            w = tgt * 2 * pi / 60
            print(f"  {tgt:6d} {got[tgt]:7.3f}  {spec.torque_at_speed(w, motor):7.3f}")
    print(f"  flat to ~{motor.corner_speed/N*60/2/pi:.0f} rpm (current-limited), "
          f"then droops to 0 at {spec.no_load_speed*60/2/pi:.0f} rpm (voltage-limited).\n")

    # --- Scenario F: drive/backdrive window holding the payload ------------------
    # The current band that keeps the loaded arm still at horizontal. Below it the
    # load backdrives the gearbox, above it the motor lifts. The legacy model puts
    # this band too high (holding costs τe/η) and too narrow (only drag resists).
    print("--- Scenario F: holding current window, payload at horizontal ---")
    model.opt.gravity[:] = (0, 0, -9.81)
    lo_s, hi_s = hold_window_sim(model, aid, fr)
    e = abs(grav_torque)
    lo_m, hi_m = (x / gear for x in fr.f.hold_window(e))
    legacy = gearbox_friction(spec, legacy=True)
    lo_l, hi_l = (x / spec.torque_per_amp for x in legacy.hold_window(e))
    print(f"  load {e:.3f} N·m   {'':14s}{'falls below':>12s} {'lifts above':>12s}")
    print(f"  this model (sim)         {lo_s*1e3:9.0f} mA {hi_s*1e3:9.0f} mA")
    print(f"  this model (analytic)    {lo_m*1e3:9.0f} mA {hi_m*1e3:9.0f} mA")
    if not LEGACY_M1:
        print(f"  legacy M1  (analytic)    {lo_l*1e3:9.0f} mA {hi_l*1e3:9.0f} mA")
        print(f"  -> minimum holding current {lo_m*1e3:.0f} mA vs {lo_l*1e3:.0f} mA legacy "
              f"({lo_l/max(lo_m,1e-9):.1f}× less): gear friction helps hold a load.")
    print()

    if "--view" in sys.argv:
        # optional initial output speed: `--view 0.05` or `--speed 0.05`
        rev = 0.10
        for flag in ("--view", "--speed"):
            if flag in sys.argv:
                i = sys.argv.index(flag)
                if i + 1 < len(sys.argv):
                    try:
                        rev = float(sys.argv[i + 1])
                    except ValueError:
                        pass
        view_demo(Params(), out_rev_per_s=rev)


def view_demo(p, out_rev_per_s=0.10):
    """Real-time, kinematically-driven look at the WHOLE hybrid (both stages).
    Each joint is driven at its true ratio so the 40:1 compound is visible:
      output 1x | cyclo carrier (= planet carrier) cyclo_ratio | sun total ratio |
      planets orbit cyclo_ratio, spin (planet_abs - carrier) about their pins.

    Speed is live-adjustable: ↑/↓ scale it, SPACE pauses, R resets. The starting
    output speed comes from `out_rev_per_s` (CLI: `--view 0.05`)."""
    import time
    from mujoco import viewer as mjv

    Nc = p.cyclo_ratio                      # output -> planet/cyclo carrier
    Nt = p.ratio                            # output -> sun (total)
    # planet absolute spin for a ring-fixed planetary (z_sun, z_planet):
    #   w_planet = w_carrier - (z_sun/z_planet)*(w_sun - w_carrier)
    wc = Nc                                 # carrier factor (rel. to output)
    ws = Nt                                 # sun factor
    w_planet_abs = wc - (p.n_sun / p.n_planet) * (ws - wc)
    planet_rel = w_planet_abs - wc          # planet spin relative to its carrier

    m = mujoco.MjModel.from_xml_path(str(HERE / "viewer_scene.xml"))
    d = mujoco.MjData(m)

    def q(name):
        return m.jnt_qposadr[m.joint(name).id]
    q_out, q_sun, q_pc = q("joint_out"), q("joint_sun"), q("joint_pcarrier")
    q_planets = [q(f"joint_planet{i}") for i in range(4)]

    # live speed control via keyboard (GLFW keycodes); mutable so the callback can edit it
    KEY_UP, KEY_DOWN, KEY_SPACE, KEY_R = 265, 264, 32, 82
    sc = {"mult": 1.0, "paused": False}

    def on_key(keycode):
        if keycode == KEY_UP:
            sc["mult"] *= 1.5
        elif keycode == KEY_DOWN:
            sc["mult"] /= 1.5
        elif keycode == KEY_SPACE:
            sc["paused"] = not sc["paused"]
        elif keycode == KEY_R:
            sc["mult"] = 1.0
        eff = out_rev_per_s * sc["mult"]
        print(f"  speed: output {eff:.3f} rev/s  (sun {eff*Nt:.2f} rev/s)"
              f"{'  [PAUSED]' if sc['paused'] else ''}")

    dt = 1.0 / 60.0
    print(f"\nviewer: output {out_rev_per_s:.2f} rev/s | carrier {out_rev_per_s*Nc:.1f} | "
          f"sun {out_rev_per_s*Nt:.1f} rev/s  ({Nt:.0f}:1 total).")
    print("  controls: ↑/↓ faster/slower · SPACE pause · R reset · close window to exit")
    with mjv.launch_passive(m, d, key_callback=on_key) as v:
        theta = 0.0
        while v.is_running():
            if not sc["paused"]:
                theta += out_rev_per_s * sc["mult"] * 2 * pi * dt
            d.qpos[q_out] = theta
            d.qpos[q_pc] = theta * Nc
            d.qpos[q_sun] = theta * Nt
            for qp in q_planets:
                d.qpos[qp] = theta * planet_rel
            mujoco.mj_forward(m, d)
            v.sync()
            time.sleep(dt)


if __name__ == "__main__":
    main()
