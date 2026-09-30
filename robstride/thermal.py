"""
Thermal model of a RobStride, fitted to the booklet's overload tables.

State is the winding's "thermal load" θ: 0 = ambient, 1 = the driver's protection
limit (145 °C, 135 °C for RS01/RS05). The rated torque is by definition what holds θ
at 1 forever, so heat input is (τ/τ_rated)² while rotating.

Holding still is worse: one phase carries the peak current continuously instead of
the RMS sharing a rotating field, 1.414× the heat per the booklet, and RobStride
publishes a lower stall rating (RS03 13 vs 21 N·m). So heat is (τ/τ_hold)² when
the rotor isn't turning. Any rotation above ~1 electrical Hz averages the phases;
that's ≈0.03 rad/s at the output, so the blend uses exp(−|ω|/0.05).

A single time constant can't fit the tables: one short burst (≈10 s at peak) and a
long tail (minutes at 1.5×) need two. Two nodes do it, winding into housing:
    θ(t) = x²·[a·(1 − e^(−t/τ1)) + (1 − a)·(1 − e^(−t/τ2))]
fitted per model; every published overload time is reproduced within ~20%
(τ1 7–26 s, τ2 3–19 min). `python thermal.py` prints the fit.
"""

import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

sys.path.insert(0, str(Path(__file__).resolve().parent))
from specs import SPECS, Spec  # noqa: E402

W_ROTATING = 0.05     # rad/s at the output: above this the phases share the heat


def _rise(t, x, a, t1, t2):
    return x * x * (a * (1 - np.exp(-t / t1)) + (1 - a) * (1 - np.exp(-t / t2)))


def time_to_limit(x, a, t1, t2):
    """Seconds from cold to θ = 1 at constant torque ratio x = τ/τ_rated."""
    if x <= 1:
        return np.inf
    return brentq(lambda t: _rise(t, x, a, t1, t2) - 1, 1e-6, 1e7)


@lru_cache(maxsize=None)
def fit(name: str):
    """(a, τ1, τ2) for SPECS[name] from its overload table."""
    s = SPECS[name]

    def res(p):
        a, l1, l2 = p
        return [np.log(time_to_limit(T / s.tau_rated, a, np.exp(l1), np.exp(l2)) / t)
                for T, t in s.overload]

    best = None
    for a0 in (0.2, 0.4, 0.6):
        r = least_squares(res, [a0, np.log(10), np.log(600)],
                          bounds=([0.02, np.log(0.5), np.log(30)],
                                  [0.98, np.log(200), np.log(5000)]))
        if best is None or r.cost < best.cost:
            best = r
    a, l1, l2 = best.x
    return float(a), float(np.exp(l1)), float(np.exp(l2))


class Thermal:
    """Two-node winding model; step(τ, ω, dt) -> θ (1 = protection trips)."""

    def __init__(self, spec: Spec):
        self.s = spec
        self.a, self.t1, self.t2 = fit(spec.name)
        self.n1 = self.n2 = 0.0

    def heat(self, tau, dq):
        w = np.exp(-abs(dq) / W_ROTATING)
        rot = (tau / self.s.tau_rated) ** 2
        hold = (tau / self.s.tau_stall_rated) ** 2
        return (1 - w) * rot + w * hold

    def step(self, tau, dq, dt):
        q = self.heat(tau, dq)
        self.n1 += dt / self.t1 * (self.a * q - self.n1)
        self.n2 += dt / self.t2 * ((1 - self.a) * q - self.n2)
        return self.theta

    @property
    def theta(self):
        return self.n1 + self.n2

    def reset(self):
        self.n1 = self.n2 = 0.0


def cycle_theta(spec: Spec, tau, dq, dt, repeats=None):
    """θ after running a periodic cycle for a long time (default: until steady, 5·τ2)
    and the peak θ reached. Faster than stepping: the slow node sees the mean heat."""
    th = Thermal(spec)
    q = np.array([th.heat(t, w) for t, w in zip(tau, dq)])
    period = len(tau) * dt
    n2 = (1 - th.a) * q.mean()          # slow node at periodic steady state
    n1, peak = th.a * q.mean(), 0.0
    for _ in range(3):                  # fast node: settle over a few cycles
        for qk in q:
            n1 += dt / th.t1 * (th.a * qk - n1)
            peak = max(peak, n1 + n2)
    return float(n1 + n2), float(peak), period


if __name__ == "__main__":
    for name, s in SPECS.items():
        a, t1, t2 = fit(name)
        pts = "  ".join(f"{T:g}N·m {t:g}s→{time_to_limit(T/s.tau_rated, a, t1, t2):.0f}s"
                        for T, t in s.overload)
        print(f"{name:6s} a={a:.2f} τ1={t1:5.1f}s τ2={t2:5.0f}s   {pts}")
