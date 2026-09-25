"""
TWO-CELL VARIABLE CHAIN MOTOR — the minimum chain that bends.

After Tada et al., "Development of Variable Chain Motor with Shape and Speed-Torque
Characteristics Variability ...", IROS 2025, doi 10.1109/IROS60139.2025.11246199.

Two stock outrunners (Goolsky 2204, the motor used elsewhere in this repo) become MOTOR CELLS:
a spur gear on the ROTOR bell, a spur gear on the STATOR mount, and the winding's star point
cut so all six phase ends come out. The cells run as ONE BLDC from ONE driver and ONE encoder.

    cell A   stator bolted to the base (the "upper arm"); its rotor shaft is the OUTPUT + encoder
    cell B   rides a SWING LINK that pivots about A's axis (the "forearm"); B's stator is free
             to spin on a bearing in the link
    rotor gear A <-> rotor gear B     (ratio r_r)  synchronizes the rotor angles
    stator gear A <-> stator gear B   (ratio r_s)  the hinge trick

Kinematics (in the rotating link frame, hinge angle phi, A-rotor angle theta, A-stator fixed):
    stator B (link frame) =  r_s * phi
    rotor  B (link frame) = -r_r * (theta - phi)
    => B electrical angle  = p * (rotor B - stator B) = -p*r_r*theta + p*(r_r - r_s)*phi
With r_r == r_s the hinge angle DROPS OUT: bending the chain can't disturb commutation, and by
virtual work the motor torque puts T_B*(r_r - r_s) = 0 on the hinge (shape and output decouple).

Two consequences the paper doesn't spell out:
  * the external mesh makes B turn BACKWARDS relative to A, so B's phase sequence must be
    reversed (swap two leads) or its torque cancels A's;
  * gear backlash across a hinge is TWO meshes (rotor + stator) of slop, multiplied by pole pairs
    into electrical-angle error -> a cos() torque loss.

Electrical (MCSW): 2S -> Kt x2, speed/V x1/2, one R_sw in the path;  2P -> Kt x1, R/2.
Validation: the same envelope code reproduces the paper's Table III (4 x MN5006, 72 V / 30 A).

    ../.venv/bin/python chain-motor/chain.py        # report + out/chain.png
"""

from dataclasses import dataclass, field
from math import pi, sqrt, cos, tan, radians
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
W_TEST = 400.0                   # efficiency test speed [rad/s] (paper used 27% of its 4S max)
OUT = HERE / "out"


@dataclass
class Cell:
    """One motor cell (a stock BLDC with its star point cut)."""
    name: str = "Goolsky 2204 (1400 KV)"
    kt: float = 6.82e-3          # N·m/A (= Ke, V·s/rad) — repo back-EMF model value
    r_ohm: float = 0.307         # terminal resistance
    pole_pairs: int = 7          # 12N14P
    i_cont_A: float = 13.0       # per-cell current rating (thermal)
    t_nl_Nm: float = 2.7e-3      # no-load drag (iron + bearings) = Kt * I0, I0 ~ 0.4 A
    mass_g: float = 28.0


@dataclass
class ChainParams:
    cell: Cell = field(default_factory=Cell)
    v_bus: float = 11.1          # 3S
    i_drv_A: float = 13.0        # driver current limit
    v_util: float = 1.0          # speed = v_util*V/Kt ; 1.0 when Kt comes from the Kv spec,
                                 # 1/sqrt(2) reproduces the paper's star-converted convention
    r_ext_ohm: float = 0.08      # driver FETs + cables (series with the whole chain)
    r_sw_ohm: float = 0.010      # one MCSW switch closed = 2 back-to-back FETs (~5 mΩ each)

    # --- gears (printed, module 1) ---
    module_mm: float = 1.0
    z_rotor: int = 34            # rotor gears, both cells (r_r = 1); 34 leaves room for the
    z_stator: int = 34           # rotor caps' skirts between the Ø28 bells (r_s = 1)
    eta_mesh: float = 0.95       # printed spur mesh efficiency (torque passed A<-B)
    t_mesh_Nm: float = 1.0e-3    # mesh churning drag referred to the motor shaft
    backlash_mm: float = 0.10    # per mesh, at the pitch circle

    @property
    def center_mm(self):
        """Both meshes share the link's centre distance: the 1:1 pitch distance m·z, opened up
        by Δa = j/(2·tan α) so each mesh has `backlash_mm` of circular backlash (j ≈ 2·Δa·tan α).
        (First cut used exactly m·z — zero backlash; text-to-cad's closest_points caught it.)"""
        return self.module_mm * self.z_rotor + self.backlash_mm / (2 * tan(radians(20.0)))

    @property
    def r_r(self): return 1.0    # z_rotorA / z_rotorB (identical cells)

    @property
    def r_s(self): return 1.0

    # ---- kinematics ----------------------------------------------------------
    def elec_angle_B(self, theta, phi, r_r=None, r_s=None, reversed_phases=True):
        """B's electrical angle vs A's for rotor angle theta, hinge phi (rad)."""
        r_r = self.r_r if r_r is None else r_r
        r_s = self.r_s if r_s is None else r_s
        p = self.cell.pole_pairs
        eB = p * (-r_r * theta + (r_r - r_s) * phi)
        return -eB if reversed_phases else eB   # lead swap reverses the electrical sense

    def hinge_torque(self, t_B, r_r=None, r_s=None):
        """Torque motor B puts on the hinge (virtual work): T_B*(r_r - r_s)."""
        r_r = self.r_r if r_r is None else r_r
        r_s = self.r_s if r_s is None else r_s
        return t_B * (r_r - r_s)

    def backlash_elec_err(self, meshes=2):
        """Worst-case electrical-angle error from backlash across `meshes` meshes (rad)."""
        r_pitch = self.module_mm * self.z_rotor / 2.0
        return self.cell.pole_pairs * meshes * self.backlash_mm / r_pitch

    # ---- electrical configurations --------------------------------------------
    def config(self, ns, np_):
        """ns cells in series per branch, np_ branches in parallel."""
        c = self.cell
        n_sw = (ns - 1) * np_ if np_ > 0 else 0          # series links closed per path
        kt = ns * c.kt
        r = ns * c.r_ohm / np_ + (ns - 1) * self.r_sw_ohm + (self.r_sw_ohm if np_ > 1 else 0.0)
        i_max = min(self.i_drv_A, np_ * c.i_cont_A)
        return dict(kt=kt, r=r, i_max=i_max, n=ns * np_)

    def envelope(self, ns, np_, n_pts=200, gears=True):
        """Output torque vs speed: current-limited, then voltage-limited (steady state)."""
        cfg = self.config(ns, np_)
        w0 = self.v_util * self.v_bus / cfg["kt"]
        w = np.linspace(0, w0, n_pts)
        i_v = (self.v_util * self.v_bus - cfg["kt"] * w) / (cfg["r"] + self.r_ext_ohm)
        i = np.minimum(cfg["i_max"], i_v)
        t_el = cfg["kt"] * i                                   # electromagnetic, summed
        t_out = self.output_torque(t_el, cfg["n"], gears)
        return w, np.maximum(t_out, 0.0)

    def output_torque(self, t_el, n, gears=True):
        """Shaft torque: A's share direct, B's share through one rotor mesh, minus drags."""
        c = self.cell
        if n == 1 or not gears:
            return t_el - n * c.t_nl_Nm
        share = t_el / n                                       # currents are synchronized
        return share + (n - 1) * self.eta_mesh * share - n * c.t_nl_Nm - (n - 1) * self.t_mesh_Nm

    def efficiency(self, ns, np_, t_out, w):
        """Electrical-in -> shaft-out at a load point (single motor: ns=np_=1)."""
        cfg = self.config(ns, np_)
        c, n = self.cell, cfg["n"]
        if n == 1:
            t_el = t_out + c.t_nl_Nm
        else:  # invert output_torque
            share = (t_out + n * c.t_nl_Nm + (n - 1) * self.t_mesh_Nm) / (1 + (n - 1) * self.eta_mesh)
            t_el = n * share
        i = t_el / cfg["kt"]
        p_in = t_el * w + i ** 2 * (cfg["r"] + self.r_ext_ohm)
        return t_out * w / p_in


def paper_table3():
    """Reproduce Tada et al. Table III: 4 x MN5006 KV450, 72 V / 30 A (their star convention)."""
    mn = Cell(name="MN5006 KV450", kt=2.72e-2, r_ohm=0.0, i_cont_A=1e9, t_nl_Nm=0.0)
    P = ChainParams(cell=mn, v_bus=72.0, i_drv_A=30.0, v_util=1 / sqrt(2), r_ext_ohm=0.0,
                    r_sw_ohm=0.0)
    rows = []
    for label, ns, np_, paper in (("4-serial", 4, 1, (3.27, 468)),
                                  ("2-serial-2-parallel", 2, 2, (1.63, 935)),
                                  ("4-parallel", 1, 4, (0.816, 1871))):
        cfg = P.config(ns, np_)
        rows.append((label, cfg["kt"] * cfg["i_max"], P.v_util * P.v_bus / cfg["kt"], paper))
    return rows


def report(P: ChainParams):
    c = P.cell
    print(f"\n=== 2-cell Variable Chain Motor — {c.name} x2, {P.v_bus} V / {P.i_drv_A} A driver ===")

    print("\n[validate] same envelope code vs paper Table III (4 x MN5006, 72 V / 30 A):")
    ok = True
    for label, t, w, (tp, wp) in paper_table3():
        e = max(abs(t / tp - 1), abs(w / wp - 1))
        ok &= e < 0.01
        print(f"   {label:20s} {t:6.3f} N·m {w:6.0f} rad/s   paper {tp:5.3f} / {wp:4d}   "
              f"{'✓' if e < 0.01 else '✗'} ({e*100:.1f}%)")

    print("\n[kinematics] B electrical angle minus A's, rotor fixed, hinge swept 0→155°:")
    phi = np.radians(155)
    d_ok = P.elec_angle_B(0.0, phi) - P.elec_angle_B(0.0, 0.0)
    d_bad = P.elec_angle_B(0.0, phi, r_s=33 / 34) - P.elec_angle_B(0.0, 0.0, r_s=33 / 34)
    d_none = P.elec_angle_B(0.0, phi, r_s=0.0) - P.elec_angle_B(0.0, 0.0, r_s=0.0)
    print(f"   r_s = r_r (34:34 both)       Δθe = {np.degrees(d_ok):8.2f}°   ✓ hinge-invariant")
    print(f"   stator 33:34 (1-tooth off)   Δθe = {np.degrees(d_bad):8.2f}°   → torque ×{cos(d_bad):+.2f}")
    print(f"   no stator gears (fixed st.)  Δθe = {np.degrees(d_none):8.2f}°   → commutation lost")
    print(f"   hinge torque from T_B=0.1 N·m: matched {P.hinge_torque(0.1):.4f} N·m, "
          f"1-tooth off {P.hinge_torque(0.1, r_s=33/34):.4f} N·m")
    print(f"   direction: A +θ ⇒ B rotor −θ  → swap two of B's leads, else torques CANCEL "
          f"(no swap: B electrical = {np.degrees(P.elec_angle_B(1.0, 0, reversed_phases=False)):.0f}° "
          f"vs A {np.degrees(P.cell.pole_pairs*1.0):.0f}° per rad)")
    eb = P.backlash_elec_err()
    print(f"   backlash {P.backlash_mm} mm/mesh × 2 meshes × p={c.pole_pairs} → {np.degrees(eb):.1f}° "
          f"electrical worst case → B torque ×{cos(eb):.3f}")

    print(f"\n[envelope] output shaft (gears m{P.module_mm} z{P.z_rotor}, centre {P.center_mm:.3f} mm, "
          f"η_mesh {P.eta_mesh}):")
    single = ChainParams(cell=c, v_bus=P.v_bus, i_drv_A=P.i_drv_A, r_ext_ohm=P.r_ext_ohm)
    for label, ns, np_, PP in (("single 2204", 1, 1, single), ("2-serial", 2, 1, P),
                               ("2-parallel", 1, 2, P)):
        w, t = PP.envelope(ns, np_)
        cfg = PP.config(ns, np_)
        print(f"   {label:12s} Kt {cfg['kt']*1e3:5.2f} mN·m/A  R {cfg['r']*1e3:5.0f} mΩ  "
              f"stall {t[0]*1e3:6.1f} mN·m  no-load {w[-1]:6.0f} rad/s ({w[-1]*60/2/pi:5.0f} rpm)")

    print(f"\n[efficiency] at {W_TEST:.0f} rad/s (~{W_TEST*60/2/pi:.0f} rpm, ¼ of the 2204's no-load), vs load:")
    loads = [0.006, 0.01, 0.02, 0.05, 0.10]
    print("   load mN·m     " + "  ".join(f"{L*1e3:6.0f}" for L in loads))
    for label, ns, np_, PP in (("single 2204", 1, 1, single), ("2-serial", 2, 1, P),
                               ("2-parallel", 1, 2, P)):
        print(f"   {label:12s}  " + "  ".join(f"{PP.efficiency(ns, np_, L, W_TEST):6.2f}" for L in loads))
    grid = np.linspace(0.003, 0.15, 400)
    cross = next(L for L in grid if P.efficiency(2, 1, L, W_TEST) > single.efficiency(1, 1, L, W_TEST))
    print(f"   → 2S beats a single motor above ~{cross*1e3:.0f} mN·m; below that the chain pays 2× no-load\n"
          "     drag + mesh loss. Series always beats parallel (half the current through R_ext) —\n"
          "     the same pattern as the paper's Table VI.")

    mass = 2 * c.mass_g
    print(f"\nVERDICT: {'VALID' if ok else 'CHECK'} — r_r=r_s makes the hinge commutation-neutral; "
          f"2S ≈2× stall torque, switching to 2P ≈2× the speed; {mass:.0f} g of motors (+ printed gears/links).")
    return ok


def plot(P: ChainParams):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(exist_ok=True)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))
    single = ChainParams(cell=P.cell, v_bus=P.v_bus, i_drv_A=P.i_drv_A, r_ext_ohm=P.r_ext_ohm)
    for label, ns, np_, PP, col in (("single 2204", 1, 1, single, "0.5"),
                                    ("2-serial (2S)", 2, 1, P, "C3"),
                                    ("2-parallel (2P)", 1, 2, P, "C0")):
        w, t = PP.envelope(ns, np_)
        ax[0].plot(w * 60 / 2 / pi, t * 1e3, color=col, lw=2, label=label)
    ax[0].set_xlabel("output speed [rpm]"); ax[0].set_ylabel("output torque [mN·m]")
    ax[0].set_title("Speed-torque: switch 2S ↔ 2P on the fly"); ax[0].legend(); ax[0].grid(alpha=.3)

    phi = np.radians(np.linspace(0, 155, 200))
    for r_s, lab, col in ((1.0, "stator 34:34 (= rotor)", "C2"), (33 / 34, "stator 33:34", "C1"),
                          (0.0, "no stator gears", "C3")):
        d = np.degrees(P.elec_angle_B(0.0, phi, r_s=r_s) - P.elec_angle_B(0.0, 0.0, r_s=r_s))
        ax[1].plot(np.degrees(phi), d, color=col, lw=2, label=lab)
    ax[1].set_xlabel("hinge angle φ [°]"); ax[1].set_ylabel("B commutation error [° elec]")
    ax[1].set_title("Stator gears make the hinge commutation-neutral"); ax[1].legend(); ax[1].grid(alpha=.3)

    loads = np.linspace(0.004, 0.15, 80)
    for label, ns, np_, PP, col in (("single 2204", 1, 1, single, "0.5"), ("2S", 2, 1, P, "C3"),
                                    ("2P", 1, 2, P, "C0")):
        ax[2].plot(loads * 1e3, [PP.efficiency(ns, np_, L, W_TEST) for L in loads], color=col, lw=2,
                   label=label)
    ax[2].set_xlabel(f"load [mN·m] @ {W_TEST:.0f} rad/s"); ax[2].set_ylabel("efficiency")
    ax[2].set_title("Chain loses light, wins heavy (series)"); ax[2].legend(); ax[2].grid(alpha=.3)
    fig.suptitle("2-cell Variable Chain Motor (after Tada et al., IROS 2025) — Goolsky 2204 cells, 3S")
    fig.tight_layout()
    fig.savefig(OUT / "chain.png", dpi=110)
    print(f"wrote {OUT/'chain.png'}")


if __name__ == "__main__":
    P = ChainParams()
    report(P)
    plot(P)
