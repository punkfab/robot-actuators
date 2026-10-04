"""
The slow-manifold recipe of Bettini, Kazemipour, Katzschmann & Haller (Nat. Commun.
2026, doi 10.1038/s41467-026-77664-0) tried on our SMA cell.

Their recipe, for a HASEL muscle joint:
  1. record a few forced trajectories under slowly varying input
  2. fit a static polynomial map  output = g(input)  (the "slow manifold")
  3. feed-forward with g^-1(reference), plus a small PI loop
It needs the actuator's own transients to be fast next to the command. They measure
that with a slowness number rho; for a first-order plant with time constant tau and a
sinusoidal reference at w rad/s, rho = w * tau. The recipe is expected to hold for
rho << 1.

Here the plant is one SMA strand from bundle.py, isometric, driven by a continuous
heating duty u in [0, 1] (fraction of the actuate power). Output is activation
(force / full pull). Two things HASELs don't have make this a harder case:
  - tau is the cooling time constant, 2.7 s: rho = 0.15 (the paper's test value)
    means a reference period of about two minutes
  - the transition is hysteretic. ThermalCell uses one ramp between Mf and Af;
    HystCell below adds the loop as a play operator on temperature. Its half-width
    is [ASSUMED] (10 C); the datasheet figures in cell.py give only Mf and Af.

Four controllers on a force reference, gains tuned per controller on a separate
calibration trajectory (as the paper does):
  PI        feedback only
  FF        g^-1(r) only
  FF+PI     the paper's controller
  lead+PI   FF plus tau * d/dt[g^-1(r)]: the first-order rate correction the paper
            lists as future work. For this plant it is the exact thermal inverse.

    python muscle/slowmanifold.py        (or: make muscle-sm)
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from bundle import ThermalCell  # noqa: E402
from cell import SMACellSpec  # noqa: E402

OUT = HERE / "out"
DT = 0.02              # SMA plant step; hasel.py's plant brings its own
LOOP_C = 10.0          # [ASSUMED] hysteresis half-width, C
POLY_ORDER = 7         # the paper's order for the joint map
RHOS = (0.15, 0.5, 1.2)
CONTROLLERS = ("PI", "FF", "FF+PI", "lead+PI")
# Feedback path [ASSUMED]: without these a first-order plant takes unlimited gain and
# PI wins trivially. Tendon force sensor noise, transport delay, and the paper's EMA.
NOISE = 0.01           # rms per sample at 50 Hz, fraction of full output
DELAY = 0.10           # s, sensing + drive update
FILTER_TAU = 0.056     # s, EMA time constant (alpha = 0.3 at 50 Hz)


class HystCell:
    """ThermalCell's thermal constants, continuous duty input, hysteretic transition.

    Play operator: T_eff trails T by `loop` on heating and leads it by `loop` on
    cooling. Activation ramps over T_eff so that heating still finishes at Af and
    cooling still finishes at Mf. loop = 0 reproduces ThermalCell.activation."""

    dt = DT

    def __init__(self, spec=None, loop=LOOP_C):
        c = ThermalCell(spec or SMACellSpec())
        s = c.spec
        self.amb, self.h, self.C, self.p = s.ambient_c, c.h, c.C, c.p_pulse
        self.tau = self.C / self.h
        self.loop = loop
        self.lo = s.martensite_finish_c + loop
        self.hi = s.austenite_finish_c - loop

    def act(self, state):
        return min(max((state[1] - self.lo) / (self.hi - self.lo), 0.0), 1.0)

    def u_range(self, pad=6.0):
        """Duty span whose steady temperature covers the transition, with padding."""
        t = np.array([self.lo - self.loop - pad, self.hi + self.loop + pad])
        return (t - self.amb) * self.h / self.p

    def start(self, a0):
        """State (T, T_eff) sitting at activation a0 on the heating branch."""
        t_eff = self.lo + a0 * (self.hi - self.lo)
        return t_eff + self.loop, t_eff

    def rest(self):
        return self.amb, self.amb

    def step(self, state, u):
        T, t_eff = state
        T += (u * self.p - self.h * (T - self.amb)) / self.C * self.dt
        return T, min(max(t_eff, T - self.loop), T + self.loop)


def smooth_random(n, rng, lo, hi, n_tones=5):
    """Sum of sinusoids with random phases, scaled to [lo, hi]. One base period = n."""
    k = np.arange(n) / n
    y = sum(rng.uniform(0.4, 1.0) / m * np.sin(2 * np.pi * m * k + rng.uniform(0, 2 * np.pi))
            for m in range(1, n_tones + 1))
    y = (y - y.min()) / (y.max() - y.min())
    return lo + y * (hi - lo)


def rho_of(y, tau, dt):
    """Paper Eq. 8 for a linear plant: rms rate of the normalized signal over |lambda|.
    The mean is removed first, so an offset doesn't flatter the number."""
    y = y - y.mean()
    return tau * np.sqrt(np.sum(np.gradient(y, dt) ** 2) / np.sum(y ** 2))


def signal(cell, rho, seed, lo, hi, periods=3):
    """A random smooth signal on [lo, hi] stretched in time to hit the target rho
    for this plant (its tau and step)."""
    n = 4000
    base = smooth_random(n, np.random.default_rng(seed), lo, hi)
    n = int(round(n * rho_of(base, cell.tau, cell.dt) / rho))
    base = smooth_random(n, np.random.default_rng(seed), lo, hi)
    return np.tile(base, periods)


def forced(cell, u):
    """Open-loop response to an input trajectory, from rest."""
    state = cell.rest()
    a = np.empty(len(u))
    for k, uk in enumerate(u.tolist()):
        state = cell.step(state, uk)
        a[k] = cell.act(state)
    return a


def learn(cell, rho=0.1, n_traj=3):
    """Steps 1–2: forced trajectories -> polynomial g, and a dense inverse table."""
    u_lo, u_hi = cell.u_range()
    us, acts = [], []
    for seed in range(n_traj):
        u = signal(cell, rho, 100 + seed, u_lo, u_hi, periods=2)
        a = forced(cell, u)
        cut = len(u) // 2                     # drop the warm-up period
        us.append(u[cut:]); acts.append(a[cut:])
    u, a = np.concatenate(us), np.concatenate(acts)
    x = lambda v: 2 * (v - u_lo) / (u_hi - u_lo) - 1
    coef = np.polyfit(x(u), a, POLY_ORDER)
    g = lambda v: np.polyval(coef, x(np.clip(v, u_lo, u_hi)))
    grid_u = np.linspace(u_lo, u_hi, 2000)
    grid_a = np.maximum.accumulate(g(grid_u))            # monotone for the inverse
    g_inv = lambda r: np.interp(r, grid_a, grid_u)
    return dict(g=g, g_inv=g_inv, u=u, a=a, grid=(grid_u, g(grid_u)))


def nmte(x, xhat):
    peak = np.max(np.abs(x))
    return np.mean(np.abs(x - xhat)) / peak if peak > 0 else float("nan")


def track(cell, model, ref, ctrl, kp=0.0, ki=0.0, delay=DELAY):
    """Closed loop on the normalized output. Returns (output, input). kp and ki are
    in units of the plant's input span (ki per second), so one grid serves any plant."""
    dt = cell.dt
    span = float(np.ptp(cell.u_range()))
    ff = model["g_inv"](ref) if ctrl != "PI" else np.zeros(len(ref))
    if ctrl == "lead+PI":
        ff = ff + cell.tau * np.gradient(ff, dt)
    use_pi = ctrl != "FF"
    state = cell.start(ref[0])
    integ, filt, nd = 0.0, float(ref[0]), int(round(delay / dt))
    alpha = 1 - np.exp(-dt / FILTER_TAU)
    # same noise density at any sample rate
    noise = np.random.default_rng(0).normal(0.0, NOISE * np.sqrt(DT / dt), len(ref)).tolist()
    ref_l, ff_l = ref.tolist(), ff.tolist()
    a, us = [0.0] * len(ref), [0.0] * len(ref)
    for k in range(len(ref)):
        a[k] = cell.act(state)
        filt += alpha * (a[max(k - nd, 0)] + noise[k] - filt)
        e = ref_l[k] - filt
        raw = ff_l[k] + (span * (kp * e + ki * integ) if use_pi else 0.0)
        u = min(max(raw, 0.0), 1.0)
        if use_pi and u == raw:               # anti-windup: integrate only when unsaturated
            integ += e * dt
        us[k] = u
        state = cell.step(state, u)
    return np.array(a), np.array(us)


def rms(ref, a):
    cut = len(ref) // 3                       # score after the first period
    return float(np.sqrt(np.mean((ref[cut:] - a[cut:]) ** 2)))


GAINS = [(kp, ki) for kp in (0.1, 0.3, 1.0, 3.0, 10.0) for ki in (0.1, 0.3, 1.0, 3.0, 10.0, 30.0, 100.0)]


def tune(cell, model, ref, ctrl, delay=DELAY):
    """Best (kp, ki) on the calibration trajectory. Error is scored on the true
    output, not the filtered measurement."""
    if ctrl == "FF":
        return 0.0, 0.0
    return min(GAINS, key=lambda g: rms(ref, track(cell, model, ref, ctrl, *g, delay=delay)[0]))


def experiment(cell, rhos=RHOS, delay=DELAY):
    """Learn the map on `cell`, then predict and track at each rho. `cell` is any
    plant with tau, dt, u_range(), rest(), start(a0), step(state, u), act(state)."""
    if not hasattr(cell, "step"):
        cell = HystCell(loop=cell)
    model = learn(cell)
    out = dict(cell=cell, model=model, pred={}, rms={}, trace={}, gains={})
    u_lo, u_hi = cell.u_range()
    for rho in rhos:
        u = signal(cell, rho, 7, u_lo, u_hi)
        a = forced(cell, u)
        cut = len(u) // 3
        out["pred"][rho] = nmte(a[cut:], model["g"](u)[cut:])
        cal = signal(cell, rho, 1, 0.1, 0.9)
        test = signal(cell, rho, 2, 0.1, 0.9)
        for ctrl in CONTROLLERS:
            gains = out["gains"][rho, ctrl] = tune(cell, model, cal, ctrl, delay)
            a, us = track(cell, model, test, ctrl, *gains, delay=delay)
            out["rms"][rho, ctrl] = rms(test, a)
            out["trace"][rho, ctrl] = (test, a, us)
    return out


def plot(results):
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    for (label, r), col in zip(results.items(), ("tab:blue", "tab:red")):
        m = r["model"]
        ax[0].plot(m["u"][::25], m["a"][::25], ".", ms=2, color=col, alpha=0.35)
        ax[0].plot(*m["grid"], color=col, lw=1.8, label=f"{label}: fitted g(u)")
    ax[0].set(xlabel="heating duty u", ylabel="activation (force / full pull)",
              title="Slow-manifold fit from forced data (ρ = 0.1)")
    ax[0].legend(fontsize=8)

    r = results["hysteresis"]
    rho = RHOS[0]
    ref = r["trace"][rho, "PI"][0]
    t = np.arange(len(ref)) * r["cell"].dt
    ax[1].plot(t, ref, "k--", lw=1, label="reference")
    for ctrl in CONTROLLERS:
        ax[1].plot(t, r["trace"][rho, ctrl][1], lw=1.1, label=ctrl)
    ax[1].set(xlabel="time (s)", ylabel="activation", xlim=(t[len(t) // 3], t[-1]),
              title=f"Tracking with hysteresis, ρ = {rho}")
    ax[1].legend(fontsize=8)

    w = 0.2
    for i, ctrl in enumerate(CONTROLLERS):
        ax[2].bar(np.arange(len(RHOS)) + (i - 1.5) * w,
                  [100 * r["rms"][rho, ctrl] for rho in RHOS], w, label=ctrl)
    ax[2].set(xticks=range(len(RHOS)), xticklabels=[f"ρ = {x}" for x in RHOS],
              ylabel="RMS error (% of full force)", title="With hysteresis: error vs slowness")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    fig.savefig(OUT / "slowmanifold.png", dpi=120)


def main():
    cell = HystCell()
    print(f"SMA strand: tau = {cell.tau:.1f} s. rho = w*tau, so a sinusoidal reference has")
    for rho in RHOS:
        print(f"  rho {rho:4.2f}  ->  period {2 * np.pi * cell.tau / rho:6.1f} s")
    results = {}
    for label, loop in (("no hysteresis", 0.0), ("hysteresis", LOOP_C)):
        r = results[label] = experiment(loop)
        print(f"\n=== {label} (loop half-width {loop:.0f} C) ===")
        print("  open-loop prediction error of g(u), NMTE: "
              + "   ".join(f"rho {rho}: {100 * r['pred'][rho]:.1f}%" for rho in RHOS))
        print(f"  tracking RMS error, % of full force\n  {'rho':>5} "
              + " ".join(f"{c:>9}" for c in CONTROLLERS))
        for rho in RHOS:
            print(f"  {rho:5.2f} " + " ".join(f"{100 * r['rms'][rho, c]:9.2f}" for c in CONTROLLERS)
                  + "    gains " + " ".join(f"{c}:{r['gains'][rho, c]}" for c in CONTROLLERS if c != "FF"))
    rho = RHOS[0]
    print(f"\n=== how much hysteresis the static map survives (rho = {rho}) ===")
    print(f"  {'loop (C)':>8} {'g(u) NMTE':>10} " + " ".join(f"{c:>9}" for c in CONTROLLERS[:3]))
    for loop in (0.0, 2.0, 4.0, 7.0, LOOP_C):
        r = results["hysteresis"] if loop == LOOP_C else experiment(loop, (rho,))
        print(f"  {loop:8.0f} {100 * r['pred'][rho]:9.1f}% "
              + " ".join(f"{100 * r['rms'][rho, c]:9.2f}" for c in CONTROLLERS[:3]))
    plot(results)
    print(f"\nwrote {OUT / 'slowmanifold.png'}")


if __name__ == "__main__":
    main()
