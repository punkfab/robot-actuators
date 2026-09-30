"""
RobStride actuator data, with where each number came from.

Sources (base https://github.com/RobStride/Product_Information/blob/main/):
  [SPEC]  灵足时代RS系列产品规格介绍(2026.09.17).pdf: spec booklet: R, L, output-side
          inertia, 48 V T-N tables, overload-time tables, drawings
  [Mxx]   Product Literature/RSxx/RSxxUser Manual260713.pdf: MIT-mode ranges
  [SYSID] github.com/T-K-233/Robstride-SysID: MuJoCo sysid of unloaded actuators
          (armature, damping, Coulomb); b and Fc trade off, so loosely identified
  [KATZ]  B. Katz, "A Low Cost Modular Actuator for Dynamic Robots", MIT MS thesis
          2018, Eq. 2.27 (dyno, ±12 rad/s): the mini-cheetah actuator (17 N·m, 6:1
          single-stage planetary, the design RobStride's MIT mode comes from):
          τ_f = 0.09 + 0.04·|τ_motor|, no speed term. Used as the STRUCTURE of the
          gearbox prior (motor-side load term, η_drive 0.96, η_back 0.96); drag per
          model from [SYSID] where measured
  [KS]    K-Scale K-Bot MJCF (github.com/kscalelabs/ksim, examples/kbot): frictionloss
          values used there, not measured
  [ASSUMED] ours, marked; replace with a bench fit (record.py -> fit)

Conventions: torques and speeds at the OUTPUT. `tn` is the 48 V (RS01: 36 V) T-N
table [SPEC], (N·m, rpm) pairs; the envelope interpolates it and scales speed with
V_bus. `overload` is (N·m, seconds to thermal protection) while rotating [SPEC],
from which thermal.py fits a two-node winding/housing model. `armature` is the
booklet's "低速端等效惯量" (equivalent inertia, low-speed side) = reflected rotor.

Not published by RobStride: gearbox efficiency, backlash, loop rates. [KATZ] gives
the measured values for the ancestor design (η ≈ 0.96, backlash ≈ 0.005 rad, the MIT
PD law run on the drive at the 40 kHz current-loop rate).
EL05 is left out (no R/inertia published).
Conflicts are noted inline; the conservative value is used.
"""

from dataclasses import dataclass, field

import numpy as np

RPM = 2 * np.pi / 60


@dataclass(frozen=True)
class Spec:
    name: str
    tau_peak: float           # N·m   [SPEC]/[Mxx]; also the MIT torque range
    tau_rated: float          # N·m   rotating, continuous [SPEC]
    tau_stall_rated: float    # N·m   continuous while HOLDING STILL (one phase carries it)
    ratio: float
    v_nom: float              # V
    v_range: tuple            # V
    i_peak: float             # A
    kt: float                 # N·m/A_rms at the output [SPEC]
    r_line: float             # Ω
    armature: float           # kg·m² reflected rotor [SPEC]
    mass: float               # kg
    od: float                 # m
    length: float             # m
    v_max: float              # rad/s MIT velocity range [Mxx]
    kp_max: float             # MIT kp range (N·m/rad)
    kd_max: float             # MIT kd range (N·m·s/rad)
    tn: tuple                 # ((N·m, rpm), ...) at v_nom [SPEC]
    overload: tuple           # ((N·m, s), ...) rotating [SPEC]
    drag: float               # N·m  Coulomb at the output
    kv: float                 # N·m·s/rad viscous
    k_load: float = 0.04      # –    motor-side load-dependent friction [KATZ]: η ≈ 0.96
    temp_limit: float = 145.0 # °C winding protection [SPEC] (135 for RS01/RS05)
    p_max: float = 4 * np.pi  # MIT position range, all models [Mxx]
    notes: str = ""

    @property
    def t_max(self):
        return self.tau_peak

    @property
    def w_noload(self):
        """No-load output speed (rad/s) at v_nom: the T-N table's zero-torque end."""
        return float(self.tn_arrays()[1][-1]) * RPM

    def tn_arrays(self):
        """(torque, rpm) ascending in speed, extended to (0 N·m, no-load speed) and
        (τ_peak, 0 rpm): below the table's lowest speed the current limit rules."""
        pts = sorted(self.tn, key=lambda p: p[1])
        tau = [self.tau_peak] + [p[0] for p in pts]
        rpm = [0.0] + [p[1] for p in pts]
        if tau[-1] > 0:                        # extrapolate the last segment to τ = 0
            (t1, n1), (t2, n2) = (tau[-2], rpm[-2]), (tau[-1], rpm[-1])
            slope = (n2 - n1) / (t2 - t1) if t2 != t1 else 0.0
            tau.append(0.0)
            rpm.append(n2 - slope * t2 if slope < 0 else n2 * 1.03)
        return np.array(tau), np.array(rpm)


SPECS = {s.name: s for s in [
    Spec("RS00", tau_peak=14, tau_rated=5, tau_stall_rated=5 * 0.84, ratio=10, v_nom=48,
         v_range=(24, 60), i_peak=15.5, kt=1.48, r_line=1.5, armature=0.001,
         mass=0.310, od=0.057, length=0.051, v_max=33, kp_max=500, kd_max=5,
         tn=((0.5, 315), (5, 270), (8, 225), (12, 150), (14, 100)),
         overload=((7, 120), (10, 18), (12, 10), (14, 5)),
         drag=0.233, kv=0.02,
         notes="SPEC summary lists a newer 17/6 N·m, 330 g revision. SYSID kv 0.054 ±129%: "
               "0.02 used. stall rating [ASSUMED] rated/√1.414 (booklet: stall heat ×1.414)."),
    Spec("RS01", tau_peak=17, tau_rated=6, tau_stall_rated=6 * 0.84, ratio=7.75, v_nom=36,
         v_range=(24, 50), i_peak=23, kt=1.22, r_line=0.55, armature=0.0042,
         mass=0.380, od=0.0785, length=0.0415, v_max=44, kp_max=500, kd_max=5,
         tn=((0.5, 305), (6, 265), (10, 228), (14, 184), (17, 135)),
         overload=((8, 840), (11, 75), (14, 20), (17, 9)),
         drag=0.159, kv=0.0035, temp_limit=135,
         notes="36 V variant of the RS02 motor; rated 7 [SPEC] vs 6 [M01]: 6 used. "
               "Single encoder (motor side). Overload/friction taken from RS02."),
    Spec("RS02", tau_peak=17, tau_rated=6, tau_stall_rated=6 * 0.84, ratio=7.75, v_nom=48,
         v_range=(24, 60), i_peak=23, kt=1.22, r_line=0.55, armature=0.0042,
         mass=0.380, od=0.0785, length=0.0415, v_max=44, kp_max=500, kd_max=5,
         tn=((0.5, 407), (7, 365), (10, 326), (14, 273), (17, 219)),
         overload=((8, 840), (11, 75), (14, 20), (17, 9)),
         drag=0.159, kv=0.0035,
         notes="rated 7 [SPEC] vs 6 [M02]: 6 used. SYSID armature 0.0137 (3× SPEC; may "
               "include output hardware). IP67 version 490 g, 81×44."),
    Spec("RS03", tau_peak=60, tau_rated=21, tau_stall_rated=13, ratio=9, v_nom=48,
         v_range=(24, 60), i_peak=43, kt=2.36, r_line=0.39, armature=0.02,
         mass=0.880, od=0.098, length=0.0541, v_max=20, kp_max=5000, kd_max=100,
         tn=((10, 188), (20, 187), (30, 178), (50, 145), (60, 120)),
         overload=((30, 393), (40, 67), (50, 34), (60, 13)),
         drag=0.2, kv=0.01,
         notes="rated 21 on a 200×200 mm Al heatsink [SPEC p21]; manual says 20. "
               "Stall rated 13 [SPEC p21]. drag [KS], kv [ASSUMED]. No-load 200 [M03]/195 [SPEC]."),
    Spec("RS04", tau_peak=120, tau_rated=35, tau_stall_rated=28.5, ratio=9, v_nom=48,
         v_range=(24, 60), i_peak=90, kt=2.1, r_line=0.16, armature=0.04,
         mass=1.420, od=0.120, length=0.0422, v_max=15, kp_max=5000, kd_max=100,
         tn=((10, 190), (40, 184), (60, 160), (80, 148), (120, 95)),
         overload=((50, 290), (70, 55), (90, 24), (110, 10.6), (120, 5)),
         drag=0.2, kv=0.01,
         notes="rated 35 / overload on 220×200 mm heatsink [SPEC p25]; 40 on 345×345. "
               "Stall rated 28.5 (345² table). drag [KS], kv [ASSUMED]."),
    Spec("RS05", tau_peak=5.5, tau_rated=1.6, tau_stall_rated=1.6 * 0.84, ratio=7.75, v_nom=48,
         v_range=(15, 60), i_peak=11, kt=0.94, r_line=2.72, armature=0.0007,
         mass=0.191, od=0.046, length=0.044, v_max=50, kp_max=500, kd_max=5,
         tn=((0.5, 477), (1.6, 450), (3, 340), (4, 275), (5.5, 70)),
         overload=((2, 438), (3, 42), (4, 13), (5.5, 3.7)),
         drag=0.309, kv=0.0, temp_limit=135,
         notes="rated 1.8 [SPEC] vs 1.6 [M05]: 1.6 used. SYSID drag 0.309 is large vs "
               "its 5.5 N·m peak; SYSID armature 0.0068 (10× SPEC)."),
    Spec("RS06", tau_peak=36, tau_rated=11, tau_stall_rated=11 * 0.84, ratio=9, v_nom=48,
         v_range=(15, 60), i_peak=57, kt=1.09, r_line=0.23, armature=0.012,
         mass=0.621, od=0.082, length=0.049, v_max=50, kp_max=5000, kd_max=100,
         tn=((5, 430), (11, 426), (20, 390), (30, 320), (36, 280)),
         overload=((17, 200), (20, 36), (25, 18), (30, 8), (36, 4)),
         drag=0.169, kv=0.0093,
         notes="LeRobot PR #4070 lists v 33 / T 23 / kp 500 from the vendor SDK; manual "
               "values used. SYSID armature 0.0165."),
    Spec("RS10P", tau_peak=42, tau_rated=14, tau_stall_rated=9.5, ratio=25, v_nom=48,
         v_range=(15, 60), i_peak=19, kt=3.73, r_line=1.41, armature=0.00625,
         mass=0.460, od=0.057, length=0.0591, v_max=13, kp_max=500, kd_max=5,
         tn=((3, 120), (8, 110), (14, 100), (25, 80), (42, 30)),
         overload=((20, 68), (25, 22), (30, 10), (42, 4)),
         drag=0.3, kv=0.02, k_load=0.08,
         notes="25:1, likely two stages: k_load 0.08 (η ≈ 0.92), drag, kv all [ASSUMED]. "
               "Kt 7.73 [M10P] vs 3.73 [SPEC]: 42 N·m / (19 A/√2) ≈ 3.1 favours 3.73."),
]}


def table():
    rows = []
    for s in SPECS.values():
        rows.append(f"{s.name:6s} {s.tau_peak:5.1f}/{s.tau_rated:4.1f} N·m  hold {s.tau_stall_rated:4.1f}  "
                    f"1:{s.ratio:<5g} {s.w_noload:5.1f} rad/s  {s.mass*1e3:5.0f} g  "
                    f"Ø{s.od*1e3:.0f}×{s.length*1e3:.0f} mm  J {s.armature:.4f}  "
                    f"kp≤{s.kp_max:g}")
    return "\n".join(rows)


if __name__ == "__main__":
    print(table())
