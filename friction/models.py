"""
Extended gearbox friction models (M1–M6) after Duclusaud, Passault, Padois & Ly,
"Extended Friction Models for the Physics Simulation of Servo Actuators", ICRA 2025
(arXiv 2410.08650), and their library github.com/Rhoban/bam (Apache-2.0).

Friction is a torque BUDGET τ_f^m. Each step the applied friction is whatever would
stop the joint, clipped to ±budget. MuJoCo's dof_frictionloss does exactly that clip,
so a model only has to compute the budget (see mujoco_friction.py).

One general form covers all six models; the named ones are parameter restrictions:

  budget = Kc + Kv|θ̇| + |Km·τm − Ke·τe|
         + s(θ̇)·[Kcs + |Kms·τm − Kes·τe| + Q]        s = exp(−|θ̇/vs|^α)  (Stribeck)
  Q      = Keq·τe²  if |τm| > |τe|,  Kmq·τm²  if |τm| < |τe|  (only when τm, τe oppose)

  M1 Coulomb-viscous      Kc, Kv
  M2 + Stribeck           Kcs, vs, α
  M3 load-dependent       Km = Ke = Kl
  M4 M2 + M3              + Kms = Kes = Kls
  M5 directional          Km ≠ Ke, Kms ≠ Kes
  M6 + quadratic          Kmq, Keq        (harmonic drives)

Everything is referred to the OUTPUT joint: τm is the ideal (lossless) motor torque
at the output, τe the external torque on the output (gravity, a load, a brake).

Drive / backdrive efficiency. In both states τm and τe oppose, so the gearbox load
is Km|τm| + Ke|τe| and (Stribeck off) the steady-state balance gives
    drive:     η_f = (1 − Km) / (1 + Ke)
    backdrive: η_b = (1 − Ke) / (1 + Km)       Ke ≥ 1  ->  self-locking
from_efficiency() inverts these, so a predicted η_f / η_b becomes an M5 model.
"""

import json
from dataclasses import dataclass, fields, replace
from pathlib import Path

import numpy as np

PARAMS = Path(__file__).resolve().parent / "params"   # BAM's pre-fitted library


@dataclass
class Friction:
    kc: float = 0.0      # N·m    Coulomb (always present)
    kv: float = 0.0      # N·m·s  viscous (-> MuJoCo dof_damping, not clipped)
    km: float = 0.0      # –      load-dependent, motor side
    ke: float = 0.0      # –      load-dependent, external side
    kcs: float = 0.0     # N·m    extra Coulomb at rest (Stribeck)
    kms: float = 0.0     # –      extra load-dependent at rest, motor side
    kes: float = 0.0     # –      extra load-dependent at rest, external side
    kmq: float = 0.0     # 1/N·m  quadratic, motor side
    keq: float = 0.0     # 1/N·m  quadratic, external side
    vs: float = 0.2      # rad/s  Stribeck velocity
    alpha: float = 1.35  # –      Stribeck curvature
    name: str = "custom"

    # ---- the model ------------------------------------------------------------------
    def budget(self, tau_m, tau_e, dq):
        """Static part of the friction budget (N·m) -> dof_frictionloss.
        The viscous Kv|θ̇| part goes to dof_damping (see mujoco_friction.py)."""
        tau_m, tau_e, dq = np.asarray(tau_m), np.asarray(tau_e), np.asarray(dq)
        out = self.kc + np.abs(self.km * tau_m - self.ke * tau_e)
        if self.kcs or self.kms or self.kes or self.kmq or self.keq:
            s = np.exp(-np.abs(dq / self.vs) ** self.alpha)
            extra = self.kcs + np.abs(self.kms * tau_m - self.kes * tau_e)
            oppose = np.sign(tau_m) != np.sign(tau_e)
            q = np.where(np.abs(tau_m) > np.abs(tau_e), self.keq * tau_e**2,
                         np.where(np.abs(tau_m) < np.abs(tau_e), self.kmq * tau_m**2, 0.0))
            out = out + s * (extra + oppose * q)
        return out

    def total(self, tau_m, tau_e, dq):
        """Full budget τ_f^m including the viscous term (for analysis/plots)."""
        return self.budget(tau_m, tau_e, dq) + self.kv * np.abs(dq)

    # ---- derived numbers ------------------------------------------------------------
    @property
    def eta_drive(self) -> float:
        """Asymptotic drive efficiency (moving, high load)."""
        return (1 - self.km) / (1 + self.ke)

    @property
    def eta_back(self) -> float:
        """Asymptotic backdrive efficiency (moving). ≤ 0 means self-locking."""
        return (1 - self.ke) / (1 + self.km)

    def eta_static(self):
        """(drive, backdrive) efficiency at breakaway (Stribeck terms fully on)."""
        km, ke = self.km + self.kms, self.ke + self.kes
        return (1 - km) / (1 + ke), (1 - ke) / (1 + km)

    def hold_window(self, tau_e: float, dq: float = 0.0):
        """Motor torques (at the output) that hold |τe| still: [min, max].
        Below min it backdrives; above max it drives. Stribeck on at dq = 0.
        Solves |τe| = τm ± budget(τm, τe) with τm, τe opposing (fixed-point)."""
        e = abs(tau_e)
        lo = hi = e
        for _ in range(60):   # budget is Lipschitz < 1 in τm for sane params -> converges
            lo = max(0.0, e - float(self.budget(lo, -e, dq)))
            hi = e + float(self.budget(hi, -e, dq))
        return lo, hi

    def load_torque(self, tau_m: float, dq: float = 1.0) -> float:
        """Largest external load (N·m) that motor torque τm drives while moving at dq:
        solves T = τm − budget(τm, −T, dq) − Kv|dq|. The 'useful torque' of the drive."""
        t = tau_m
        for _ in range(60):
            t = tau_m - float(self.total(tau_m, -t, dq))
        return max(0.0, t)

    # ---- constructors ---------------------------------------------------------------
    @classmethod
    def m1(cls, kc, kv=0.0):
        return cls(kc=kc, kv=kv, name="M1")

    @classmethod
    def from_efficiency(cls, eta_drive, drag, kv=0.0, eta_back=None, **stribeck):
        """M5 prior from a predicted efficiency.
        eta_drive : asymptotic drive efficiency η∞ (e.g. cycloidal/efficiency.py)
        drag      : no-load drag at the output, N·m (what the old M1 used as Kc)
        eta_back  : backdrive efficiency; default = eta_drive (symmetric, = M3).
                    BAM's fitted servos span η_b/η_f ≈ 0.94–1.17, so symmetric is the
                    unbiased prior until the split is measured.
        Kc is set to drag·(1+Ke) so the DRIVE curve is unchanged from the old
        gear=η∞ + frictionloss=drag model: η(T) = η∞·T / (T + drag)."""
        ef = eta_drive
        eb = ef if eta_back is None else eta_back
        det = 1 - ef * eb
        km = (1 - 2 * ef + ef * eb) / det
        ke = (1 - 2 * eb + ef * eb) / det
        if km < -1e-9 or ke < -1e-9:
            raise ValueError(f"η_f={ef}, η_b={eb} need a negative coefficient "
                             f"(Km={km:.3f}, Ke={ke:.3f}); not a passive gearbox")
        km, ke = max(km, 0.0), max(ke, 0.0)
        return cls(kc=drag * (1 + ke), kv=kv, km=km, ke=ke, name="M5-prior", **stribeck)

    @classmethod
    def from_bam(cls, actuator: str, model: str = "m6"):
        """Load a BAM pre-fitted model, e.g. from_bam('sts3215', 'm4').
        Returns (Friction, extras) where extras holds kt, R, armature, etc."""
        d = json.loads((PARAMS / actuator / f"{model}.json").read_text())
        kl = d.get("load_friction_base", 0.0)          # M3/M4: one symmetric coefficient
        kls = d.get("load_friction_stribeck", 0.0)
        f = cls(
            kc=d.get("friction_base", 0.0),
            kv=d.get("friction_viscous", 0.0),
            km=d.get("load_friction_motor", kl),
            ke=d.get("load_friction_external", kl),
            kcs=d.get("friction_stribeck", 0.0),
            kms=d.get("load_friction_motor_stribeck", kls),
            kes=d.get("load_friction_external_stribeck", kls),
            kmq=d.get("load_friction_motor_quad", 0.0),
            keq=d.get("load_friction_external_quad", 0.0),
            vs=d.get("dtheta_stribeck", 0.2),
            alpha=d.get("alpha", 1.35),
            name=f"{actuator}/{model}",
        )
        extras = {k: v for k, v in d.items() if k not in {fl.name for fl in fields(cls)}}
        return f, extras

    def scaled(self, **kw):
        return replace(self, **kw)

    def describe(self) -> str:
        es, bs = self.eta_static()
        lock = "  SELF-LOCKING" if self.eta_back <= 0 else ""
        return (f"{self.name}: Kc={self.kc*1e3:.1f} mN·m Kv={self.kv:.4f} Km={self.km:.3f} "
                f"Ke={self.ke:.3f} | η drive {self.eta_drive:.0%} back {self.eta_back:.0%}"
                f" (static {es:.0%}/{bs:.0%}){lock}")


def library_table():
    """η_drive / η_back of every BAM-fitted M5 — the empirical spread for priors."""
    rows = []
    for d in sorted(p for p in PARAMS.iterdir() if (p / "m5.json").exists()):
        f, _ = Friction.from_bam(d.name, "m5")
        es, bs = f.eta_static()
        rows.append((d.name, f.km, f.ke, f.eta_drive, f.eta_back, es, bs))
    return rows


if __name__ == "__main__":
    print("BAM pre-fitted M5 models (github.com/Rhoban/bam), output-referred:")
    print(f"  {'actuator':18s} {'Km':>6s} {'Ke':>6s} {'η_f':>5s} {'η_b':>5s} {'η_b/η_f':>7s}"
          f"  {'static f/b':>10s}")
    for name, km, ke, ef, eb, es, bs in library_table():
        print(f"  {name:18s} {km:6.3f} {ke:6.3f} {ef:5.0%} {eb:5.0%} {eb/ef:7.2f}  {es:4.0%}/{bs:4.0%}")
