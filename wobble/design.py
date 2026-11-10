"""
Nutating ("wobble") roller-cam reducer whose cam rollers are also the output bearing.

The mechanism (Bostan's precessional drive, the family TDI's coaxial roller-cam reducer belongs to):
a block nutates about a fixed point O, tilted by THETA, the tilt direction turning with the input
shaft (angle phi). It carries two rows of rollers on radial pins whose axes pass through O. Each
row runs in a wavy groove cut in a ring around the block: the FIXED ring (z1 lobes, block row z2
rollers) and the OUTPUT ring (z4 lobes, block row z3 rollers). Lobe counts differ from roller
counts by one, so

    block spin   psi   = phi (1 - z1/z2)
    output       psi_o = phi (1 - z1 z3 / (z2 z4))       i = z2 z4 / (z2 z4 - z1 z3)

15/14 and 13/14 give 196:1, the ratio in TDI's demo.

The crossed-roller idea. A roller in a closed groove would touch both walls and skid. As in a
crossed roller bearing, alternate rollers are made to touch alternate walls: A rollers sit EPS
below the row's cone and ride the lower wall, B rollers sit EPS above and ride the upper wall.
Each roller only ever rolls on one surface, and the row as a whole pushes both ways. Preloaded
by a small interference, each row is a four-point-style roller bearing whose race is a cam: with
the rings split at the equator, the two rings are the halves of a crossed roller bearing's outer
ring and the block is its inner ring. The question this file answers is whether that bearing
carries the actuator's output loads (thrust both ways, tilting moment, radial) on top of the
drive torque, without a separate main bearing.

Everything is on spheres about O: rollers are cones with apex at O, so the analysis is done on
the roller row's pitch sphere (radius R, the middle of the roller) and scales along the roller.

    python wobble/design.py           (or: make wobble)   -> out/*.png, design.json
"""

import json
from dataclasses import asdict, dataclass, field, replace
from math import gcd
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"

E_STAR = 210e3 / (2 * (1 - 0.3 ** 2))          # MPa, steel on steel
P_DYN = 1500.0          # MPa  contact pressure for long life, case-hardened (gear practice, ~1e9 cycles)
P_STAT = 4000.0         # MPa  static limit for roller contacts (ISO 76)


@dataclass
class Row:
    rollers: int            # rollers on the block in this row (even: A and B alternate)
    lobes: int              # lobes on the ring's groove (rollers +- 1)
    elev: float             # deg, the row's cone: roller axes at this elevation in the block frame


@dataclass
class Wobble:
    R: float = 40.0                 # mm  pitch sphere radius (middle of the roller)
    L: float = 12.0                 # mm  roller length along its axis
    r: float = 2.75                 # mm  roller radius at R (cone, scales with radius)
    theta: float = 2.5              # deg nutation angle
    eps: float = 1.2                # deg A/B elevation split (A at elev-eps, B at elev+eps)
    preload_N: float = 150.0        # N   roller force from the preload interference, no external load
    h_aux: float = 25.0             # mm  auxiliary radial bearing on the output, this far above O
    k_aux: float = 0.0              # N/mm its radial stiffness (0: none; the twin doesn't need one)
    # (z of centre mm, clocking deg, mirrored) per block. The default twin: two blocks tilting
    # opposite ways, mirrored so both output rows sit in the middle (fixed | out | out | fixed)
    blocks: tuple = ((-20.0, 0.0, False), (20.0, 0.0, True))
    fixed: Row = field(default_factory=lambda: Row(10, 11, -16.0))
    output: Row = field(default_factory=lambda: Row(12, 13, 16.0))

    @property
    def preload(self):
        """interference (mm) per roller that gives preload_N"""
        return (self.preload_N / k_contact(self)) ** 0.9

    @property
    def ratio(self):
        z1, z2, z3, z4 = self.fixed.lobes, self.fixed.rollers, self.output.rollers, self.output.lobes
        return z2 * z4 / (z2 * z4 - z1 * z3)

    @property
    def gamma(self):
        return np.arcsin(self.r / self.R)   # roller cone half-angle


# ---------------------------------------------------------------- kinematics

def Rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def Rx(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def block_rot(w: Wobble, phi):
    """block orientation: tilt THETA about the horizontal line at azimuth phi, then spin psi"""
    psi = phi * (1 - w.fixed.lobes / w.fixed.rollers)
    return Rz(phi) @ Rx(np.radians(w.theta)) @ Rz(-phi) @ Rz(psi)


def output_angle(w: Wobble, phi):
    return phi / w.ratio


def axes(w: Wobble, row: Row):
    """roller axis unit vectors in the block frame, and A/B flags (+1 rides the lower wall)"""
    n = row.rollers
    beta = 2 * np.pi * np.arange(n) / n
    side = np.where(np.arange(n) % 2 == 0, 1, -1)                 # +1 A, -1 B
    el = np.radians(row.elev) - side * np.radians(w.eps)          # A below the cone, B above
    u = np.stack([np.cos(beta) * np.cos(el), np.sin(beta) * np.cos(el), np.sin(el)], axis=1)
    return u, side


def ring_rot(w: Wobble, row: Row, phi):
    return np.eye(3) if row is w.fixed else Rz(output_angle(w, phi))


def path(w: Wobble, row: Row, phi, u):
    """roller axis direction in its ring's frame at input angle(s) phi"""
    phi = np.atleast_1d(phi)
    return np.array([ring_rot(w, row, p).T @ block_rot(w, p) @ u for p in phi])


def frame(w: Wobble, row: Row, phi, u, h=1e-5):
    """at one phi: axis p, path tangent t (in the ring frame), and n_up (the sphere-tangent
    normal to the path pointing to higher elevation), plus geodesic curvature k (1/rad, + when
    the path bends toward n_up) -- all in the RING frame"""
    P = path(w, row, [phi - h, phi, phi + h], u)
    p = P[1]
    d1 = (P[2] - P[0]) / (2 * h)
    d2 = (P[2] - 2 * P[1] + P[0]) / h ** 2
    t = d1 / np.linalg.norm(d1)
    n = np.cross(p, t)
    up = np.array([0, 0, 1.0]) - p * p[2]
    if n @ up < 0:
        n = -n
    k = (d2 @ n) / (d1 @ d1)
    return p, t, n, k


# ---------------------------------------------------------------- groove geometry

def groove(w: Wobble, row: Row, samples=4000):
    """the two walls of a row's groove on the pitch sphere, as (azimuth, elevation) in degrees,
    plus the tightest wall curvature radius met by each roller family (mm, at R). A walls are
    the lower wall of the A axis path, B the upper wall of the B path."""
    u, side = axes(w, row)
    out = {}
    g = w.gamma
    span = 2 * np.pi * row.rollers / gcd(row.rollers, row.lobes)   # roller 0's path closes after this
    phis = np.linspace(0, span, samples, endpoint=False)
    for name, j in (("A", 0), ("B", 1)):
        pts, rad = [], []
        for ph in phis:
            p, t, n, k = frame(w, row, ph, u[j])
            s = side[j]
            wall = np.cos(g) * p - s * np.sin(g) * n
            pts.append(wall)
            # roller on its wall: path bending toward the roller's wall side (k*s < 0) makes the
            # wall convex there, radius (1/|k| - gamma); the other way concave (1/|k| + gamma)
            kk = k * s
            rad.append(w.R * ((1 / abs(kk) - g) if kk < 0 else (1 / max(abs(kk), 1e-9) + g)) * (-1 if kk >= 0 else 1))
        pts = np.array(pts)
        az = np.degrees(np.unwrap(np.arctan2(pts[:, 1], pts[:, 0])))
        el = np.degrees(np.arcsin(np.clip(pts[:, 2], -1, 1)))
        rad = np.array(rad)                       # + convex (mm), - concave
        out[name] = dict(az=az, el=el, convex_min=float(rad[rad > 0].min()) if (rad > 0).any() else np.inf)
    return out


def undercut_limit(w: Wobble, row: Row):
    """largest roller radius at R before the convex wall cusps (path curvature radius at R)"""
    u, side = axes(w, row)
    span = 2 * np.pi * row.rollers / gcd(row.rollers, row.lobes)
    ks = [abs(frame(w, row, ph, u[0])[3]) for ph in np.linspace(0, span, 3000, endpoint=False)]
    return w.R / max(ks)


def flank_angle(w: Wobble, row: Row):
    """steepest path slope (deg): the pressure angle range a roller sweeps"""
    u, _ = axes(w, row)
    span = 2 * np.pi * row.rollers / gcd(row.rollers, row.lobes)
    a = []
    for ph in np.linspace(0, span, 2000, endpoint=False):
        p, t, n, k = frame(w, row, ph, u[0])
        az = np.array([-p[1], p[0], 0]) / np.hypot(p[0], p[1])
        a.append(np.degrees(np.arctan2(abs(t[2]) / np.sqrt(1 - p[2] ** 2), abs(t @ az))))
    return max(a)


# ---------------------------------------------------------------- contact statics

def k_contact(w: Wobble):
    """roller force from approach, Palmgren line contact, two contacts in series (wall + pin)"""
    return (w.L ** 0.8 / (2 * 3.84e-5)) ** (1 / 0.9)          # Q = K * delta^(10/9), N and mm


def placements(w: Wobble):
    """each block as (A, o): the canonical block (centre O, fixed row below, output row above)
    mapped by the linear map A (a clocking about z, optionally a mirror z -> -z) and moved to o.
    A mirrored block tilts the opposite way at the same input angle and keeps the same ratio."""
    out = []
    for z0, clock, mirror in w.blocks:
        A = Rz(np.radians(clock)) @ np.diag([1, 1, -1.0 if mirror else 1.0])
        out.append((A, np.array([0, 0, z0])))
    return out


def contacts(w: Wobble, phi):
    """every roller contact at input angle phi, in WORLD coordinates (mm): block index, row,
    A/B side, point C, normal n (into the roller), curvature k, Jacobian J (delta = preload + J x),
    x = (v, rot) of each block, then (v, rot) of the output; rotations about the world origin"""
    out = []
    g = w.gamma
    nb = len(w.blocks)
    for bi, (A, o) in enumerate(placements(w)):
        for row in (w.fixed, w.output):
            u, side = axes(w, row)
            Rr = ring_rot(w, row, phi)
            for j in range(row.rollers):
                p, t, n, k = frame(w, row, phi, u[j])
                p, n = A @ Rr @ p, A @ Rr @ n                          # to world
                s = side[j]
                C = w.R * (np.cos(g) * p - s * np.sin(g) * n) + o       # contact point on the wall
                # roller in block, wall in ring; d_rel = roller - wall; A: delta = pre - d_rel.n, B: + d_rel.n
                Jb = -s * np.concatenate([n, np.cross(C, n)])
                J = np.zeros(6 * (nb + 1))
                J[6 * bi:6 * bi + 6] = Jb
                if row is w.output:
                    J[6 * nb:] = -Jb
                out.append(dict(block=bi, row="fixed" if row is w.fixed else "output", side=int(s), C=C,
                                n=n * s, k=k * s, J=J))
    return out


def contact_arrays(w: Wobble, phi, h=1e-5):
    """contacts() vectorised for the simulation: arrays over every roller of every block.
    C (N,3) contact point, mm, world; n (N,3) unit normal the wall pushes the roller along;
    t (N,3) unit path tangent (roller motion relative to its wall as phi increases);
    v (N,) roller travel along its wall per radian of input (mm/rad);
    k (N,) path curvature toward the roller's wall (for hertz); block (N,) index; out (N,) bool
    output row. Rings at their ideal pose for phi."""
    g = w.gamma
    Cs, ns, ts, ks, bs, os_ = [], [], [], [], [], []
    rows = []
    for row in (w.fixed, w.output):
        u, side = axes(w, row)
        P = np.stack([(ring_rot(w, row, p).T @ block_rot(w, p)) @ u.T for p in (phi - h, phi, phi + h)])  # 3,3,n
        p = P[1]
        d1 = (P[2] - P[0]) / (2 * h)
        d2 = (P[2] - 2 * P[1] + P[0]) / h ** 2
        sp = np.sqrt((d1 * d1).sum(0))
        t = d1 / sp
        n = np.cross(p.T, t.T).T
        up = np.array([[0], [0], [1.0]]) - p * p[2]
        flip = np.where((n * up).sum(0) < 0, -1.0, 1.0)
        n = n * flip
        k = (d2 * n).sum(0) / sp ** 2
        Rr = ring_rot(w, row, phi)
        rows.append((row is w.output, Rr @ p, Rr @ t, Rr @ n, k, side, w.R * sp))
    vs = []
    for bi, (A, o) in enumerate(placements(w)):
        for is_out, p, t, n, k, side, v in rows:
            vs.append(v)
            pw, tw, nw = A @ p, A @ t, A @ n
            C = w.R * (np.cos(g) * pw - side * np.sin(g) * nw) + o[:, None]
            Cs.append(C.T)
            ns.append((nw * side).T)
            ts.append(tw.T)
            ks.append(k * side)
            bs.append(np.full(len(k), bi))
            os_.append(np.full(len(k), is_out))
    return dict(C=np.concatenate(Cs), n=np.concatenate(ns), t=np.concatenate(ts), k=np.concatenate(ks),
                block=np.concatenate(bs), out=np.concatenate(os_), v=np.concatenate(vs))


def hertz(w: Wobble, Q, k):
    """peak line-contact pressure (MPa). k: path curvature toward the roller's own wall (1/rad);
    +k makes the wall concave (conforming), -k convex"""
    R, g = w.R, w.gamma
    r = w.r
    wall = (1 / np.maximum(np.abs(k), 1e-9) + np.sign(k) * g) * R * np.sign(k)   # + concave radius
    inv = 1 / r - 1 / wall                                       # concave subtracts, convex adds
    inv = np.maximum(inv, 1e-6)
    return np.sqrt(np.maximum(Q, 0) * E_STAR * inv / (np.pi * w.L))


def block_spin_axis(w: Wobble, phi, h=1e-6):
    """canonical block angular velocity per unit input angle (rad/rad)"""
    A = block_rot(w, phi + h) @ block_rot(w, phi - h).T
    W = (A - A.T) / 2 / (2 * h)
    return np.array([W[2, 1], W[0, 2], W[1, 0]])


def solve(w: Wobble, phi=0.0, T=0.0, Fa=0.0, Fr=0.0, M=0.0, h_load=30.0, k_in=1e10, h_aux=None, k_aux=None):
    """static equilibrium with the input locked. Loads on the output: torque T (N mm, about z),
    axial Fa (N, +z), radial Fr (N, along x) acting h_load mm above the origin, tilting moment M
    (N mm, about x). The input shaft holds only each block's tilt direction (rotations of the
    block normal to its own axis); blocks and output are otherwise located by the roller contacts,
    plus an optional auxiliary radial bearing on the output h_aux above the origin.

    Each roller row leaves one motion free: its own conjugate rolling. For a fixed row that is
    the drive, held by the input. For an output row it is a coupled tilt-and-spin wobble of the
    output about the block's centre, which the input does not hold. One block: only a support
    away from the centre stops it (the aux bearing). Two blocks tilted opposite: the two output
    wobbles disagree, so the rollers lock the output by themselves."""
    h_aux = w.h_aux if h_aux is None else h_aux
    k_aux = w.k_aux if k_aux is None else k_aux
    cs = contacts(w, phi)
    nb = len(w.blocks)
    nd = 6 * (nb + 1)
    io = 6 * nb                                              # output dof offset
    J = np.array([c["J"] for c in cs])
    K = k_contact(w)
    pl = placements(w)
    a_b = block_rot(w, phi) @ np.array([0, 0, 1.0])
    Pin = np.zeros((nd, nd))                                 # input constraint on block tilt rotations
    for bi, (A, o) in enumerate(pl):
        ab = A @ a_b
        Pin[6 * bi + 3:6 * bi + 6, 6 * bi + 3:6 * bi + 6] = np.eye(3) - np.outer(ab, ab)
    G = np.zeros(nd)
    G[io:io + 3] = [Fr, 0, Fa]
    G[io + 3:] = np.cross([0, 0, h_load], [Fr, 0, 0]) + np.array([M, 0, T])
    sx = np.tile(np.concatenate([np.ones(3), np.full(3, w.R)]), nb + 1)   # rotations as mm at R
    Ja = np.zeros((2, nd))                                   # aux bearing: x, y at (0, 0, h_aux)
    Ja[0, io], Ja[1, io + 1] = 1, 1
    Ja[0, io + 4], Ja[1, io + 3] = h_aux, -h_aux             # r x (0, 0, h) = (ry h, -rx h, 0)

    def f(y):
        x = y / sx
        dp = np.maximum(w.preload + J @ x, 0)
        Q = K * dp ** (10 / 9)
        ua = Ja @ x
        U = (K * 9 / 19 * dp ** (19 / 9)).sum() - G @ x + 0.5 * k_in * x @ Pin @ x + 0.5 * k_aux * ua @ ua
        g = J.T @ Q - G + k_aux * Ja.T @ ua + k_in * Pin @ x
        return U, g / sx

    def hess(y):
        x = y / sx
        d = np.maximum(w.preload + J @ x, 0)
        H = (J * (K * 10 / 9 * d ** (1 / 9))[:, None]).T @ J + k_in * Pin + k_aux * Ja.T @ Ja
        return H / np.outer(sx, sx)

    y = np.zeros(nd)
    U, g = f(y)
    for _ in range(200):
        dy = -np.linalg.solve(hess(y) + 1e-9 * np.eye(nd), g)
        step = 1.0
        while True:
            U2, g2 = f(y + step * dy)
            if U2 <= U + 1e-4 * step * (g @ dy) or step < 1e-8:
                break
            step /= 2
        y, U, g = y + step * dy, U2, g2
        if np.abs(g).max() < 1e-6 * (1 + np.abs(G).max()):
            break
    x = y / sx
    Q = K * np.maximum(w.preload + J @ x, 0) ** (10 / 9)
    p = hertz(w, Q, np.array([c["k"] for c in cs]))
    Om = block_spin_axis(w, phi)
    tau_in = 0.0
    for bi, (A, o) in enumerate(pl):
        M_in = -k_in * Pin[6 * bi + 3:6 * bi + 6, 6 * bi + 3:6 * bi + 6] @ x[6 * bi + 3:6 * bi + 6]
        tau_in += M_in @ (np.linalg.det(A) * A @ Om)         # angular velocity is a pseudovector
    return dict(x=x, xo=x[io:], Q=Q, p=p, F_aux=-k_aux * Ja @ x, rows=[c["row"] for c in cs],
                block=[c["block"] for c in cs], side=[c["side"] for c in cs], tau_in=tau_in,
                ok=bool(np.abs(g).max() < 1e-3 * (1 + np.abs(G).max())))


def worst(w: Wobble, phases=12, **load):
    """the load case over a roller pitch of input angle: peak pressure, peak and min roller load"""
    pitch = 2 * np.pi * w.fixed.rollers / w.fixed.rollers      # one input turn moves every roller one lobe
    out = None
    for ph in np.linspace(0, pitch, phases, endpoint=False):
        s = solve(w, ph, **load)
        if out is None or s["p"].max() > out["p"].max():
            out = s
        out.setdefault("unloaded", 0)
        out["unloaded"] = max(out["unloaded"], int((s["Q"] <= 0).sum()))
    return out


def capacity(w: Wobble, key, p_lim=P_DYN, hi=1e6, unit=1.0, **base):
    """largest load (one component) before peak pressure hits p_lim, or first roller unloads
    if 'unload' is passed as p_lim"""
    def bad(v):
        s = worst(w, phases=6, **{**base, key: v * unit})
        return (s["unloaded"] > 0) if p_lim == "unload" else (s["p"].max() > p_lim)
    lo, up = 0.0, 1.0
    while not bad(up):
        lo, up = up, up * 2
        if up > hi:
            return np.inf
    for _ in range(18):
        mid = (lo + up) / 2
        lo, up = (mid, up) if not bad(mid) else (lo, mid)
    return lo


def stiffness(w: Wobble, load=50e3):
    """output tilt and torsional stiffness (N m / rad) and radial (N / um), about preload"""
    a = solve(w, 0.0, M=load)["xo"][3]
    b = solve(w, 0.0, T=load)["xo"][5]
    c = solve(w, 0.0, Fr=1000.0, h_load=0.0)["xo"][0]
    return (load / 1e3) / abs(a), (load / 1e3) / abs(b), 1000.0 / abs(c) / 1e3


def single(w: Wobble):
    """the one-block version with its auxiliary output bearing, for comparison"""
    return replace(w, blocks=((0.0, 0.0, False),), k_aux=1e5)


# reference: the crossed roller main bearing of a Harmonic Drive HPF-25 (80 mm pitch circle)
HPF25 = dict(name="HPF-25 crossed roller", dp=80, C=11400, M_allow=410, K_tilt=37.9e4,
             src="https://www.harmonicdrive.net/products/servo-mount-gearheads/harmonic-planetary-gear-units/hpf/25/hpf-25a-11-f0u1")


def efficiency(w: Wobble, T, mu, phases=8):
    """forward efficiency at output torque T (N mm): per input radian the output takes T/|i| and
    each roller loses mu * Q * (its travel along the wall). Rolling friction mu is the tangential
    drag per unit roller load (stem friction scaled by stem/roller radius, plus rolling
    resistance). Preload drag does not shrink with load, so efficiency is load-dependent like a
    preloaded bearing's."""
    out = []
    for ph in np.linspace(0, 2 * np.pi, phases, endpoint=False):
        s = solve(w, ph, T=T)
        v = contact_arrays(w, ph)["v"]
        loss = mu * (s["Q"] * v).sum()
        work = abs(T) / abs(w.ratio)
        out.append(work / (work + loss))
    return float(np.mean(out))


# ---------------------------------------------------------------- ratio options

def ratio_table(lo=40, hi=220, zmin=10, zmax=22):
    rows = []
    for z2 in range(zmin, zmax + 1, 2):
        for z1 in (z2 - 1, z2 + 1):
            for z3 in range(zmin, zmax + 1, 2):
                for z4 in (z3 - 1, z3 + 1):
                    den = z2 * z4 - z1 * z3
                    if den == 0:
                        continue
                    i = z2 * z4 / den
                    if lo <= abs(i) <= hi:
                        rows.append((abs(i), z1, z2, z3, z4))
    rows.sort()
    return rows


# ---------------------------------------------------------------- main

def report(w: Wobble, name):
    caps = {}
    for key, unit, label in (("T", 1e3, "torque_Nm"), ("Fa", 1.0, "axial_N"), ("Fa", -1.0, "axial_reverse_N"),
                             ("M", 1e3, "moment_Nm"), ("Fr", 1.0, "radial_N")):
        caps[label] = capacity(w, key, unit=unit, h_load=0.0)
    kt, kq, kr = stiffness(w)
    s = solve(w, 0.3, T=100e3)
    out = dict(capacity_at_1500MPa=caps, tilt_Nm_rad=kt, torsion_Nm_rad=kq, radial_N_um=kr,
               tau_in_check=(s["tau_in"] / 1e3, 100 / w.ratio), aux_N_at_100Nm=float(np.abs(s["F_aux"]).max()))
    print(f"-- {name}")
    print("   capacity at %.0f MPa (worst input phase): " % P_DYN +
          ", ".join(f"{k} {v:.0f}" for k, v in caps.items()))
    print(f"   stiffness: tilt {kt:.3g} N m/rad, torsion {kq:.3g} N m/rad, radial {kr:.0f} N/um")
    print(f"   torque check (virtual work): 100 N m out needs {abs(s['tau_in']) / 1e3:.4f} N m in, expected 100/|i| = {100 / abs(w.ratio):.4f}"
          + (f"; aux bearing {out['aux_N_at_100Nm']:.0f} N" if w.k_aux else ""))
    return out


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(exist_ok=True)
    w = Wobble()
    print(f"ratio {w.ratio:.1f}:1  (fixed {w.fixed.lobes} lobes / {w.fixed.rollers} rollers,"
          f" output {w.output.lobes} lobes / {w.output.rollers} rollers), theta {w.theta} deg, R {w.R} mm")

    sym = []
    for row in (w.fixed, w.output):
        u, _ = axes(w, row)
        P0 = path(w, row, np.linspace(0, 2 * np.pi * row.rollers, 40000), u[0])
        Pj = path(w, row, np.linspace(0, 2 * np.pi, 40), u[2])
        sym.append(max(np.min(np.linalg.norm(P0 - q, axis=1)) for q in Pj) * w.R)
    print(f"groove check: roller 2 lies on roller 0's groove within {max(sym) * 1e3:.1f} um (sampling)")
    lim = {nm: undercut_limit(w, row) for nm, row in (("fixed", w.fixed), ("output", w.output))}
    fl = {nm: flank_angle(w, row) for nm, row in (("fixed", w.fixed), ("output", w.output))}
    print(f"roller radius {w.r} mm; undercut limit fixed {lim['fixed']:.2f}, output {lim['output']:.2f} mm;"
          f" steepest flank fixed {fl['fixed']:.1f}, output {fl['output']:.1f} deg")
    print(f"preload {w.preload_N:.0f} N per roller = {w.preload * 1e3:.2f} um interference")
    one = single(w)
    print(f"one block: the output's free wobble needs a tilt reaction of T / (z4 sin theta) = "
          f"{1 / (w.output.lobes * np.sin(np.radians(w.theta))):.2f} x output torque")
    r1 = report(one, "one block + aux radial bearing 25 mm above the centre")
    r2 = report(w, "twin, opposite tilt, no bearing but the rollers")
    print("   efficiency (forward, steady): preload x roller friction x output torque")
    eff = {}
    for pre in (50.0, 150.0):
        for mu, mname in ((0.004, "needle"), (0.016, "plain stem")):
            ww = replace(w, preload_N=pre)
            row = [efficiency(ww, T * 1e3, mu) for T in (10, 40, 80)]
            eff[f"preload {pre:.0f} N, {mname} mu {mu}"] = row
            print(f"     preload {pre:3.0f} N, {mname:10s} (mu {mu}): " +
                  ", ".join(f"{T} N m {e * 100:4.1f}%" for T, e in zip((10, 40, 80), row)))
    r2["efficiency"] = eff
    print(f"-- reference {HPF25['name']}: C {HPF25['C']} N, allowable moment {HPF25['M_allow']} N m, "
          f"tilt {HPF25['K_tilt']:.3g} N m/rad")

    # figures: the grooves on the pitch sphere, and roller loads under a few cases
    fig, ax = plt.subplots(2, 1, figsize=(11, 6.5), sharex=True)
    for a, (nm, row) in zip(ax, (("fixed", w.fixed), ("output", w.output))):
        gv = groove(w, row)
        u, side = axes(w, row)
        for key, col in (("A", "#d9480f"), ("B", "#1c7ed6")):
            az = np.mod(gv[key]["az"], 360)
            o = np.argsort(az)
            a.plot(az[o], gv[key]["el"][o], ".", ms=1, color=col,
                   label=f"{key} wall ({'lower' if key == 'A' else 'upper'})")
        P = np.array([frame(w, row, 0.0, uj)[0] for uj in u])
        a.scatter(np.mod(np.degrees(np.arctan2(P[:, 1], P[:, 0])), 360), np.degrees(np.arcsin(P[:, 2])),
                  s=50, c=["#d9480f" if sd > 0 else "#1c7ed6" for sd in side], edgecolors="k", zorder=3)
        a.set_ylabel("elevation (deg)")
        a.set_title(f"{nm} ring: {row.lobes}-lobe groove, {row.rollers} rollers (dots = roller axes at phi = 0)")
        a.legend(loc="upper right", fontsize=8)
    ax[1].set_xlabel("azimuth (deg)")
    fig.tight_layout()
    fig.savefig(OUT / "grooves.png", dpi=130)

    cases = [("preload", {}), ("T 100 N m", dict(T=100e3)), ("axial 3 kN", dict(Fa=3e3)),
             ("moment 150 N m", dict(M=150e3))]
    fig, ax = plt.subplots(1, 2, figsize=(12, 4))
    for nm, ld in cases:
        s = solve(w, 0.0, **ld)
        rows, blk = np.array(s["rows"]), np.array(s["block"])
        for a, rn in zip(ax, ("fixed", "output")):
            m = rows == rn
            a.plot(np.arange(m.sum()), s["Q"][m], "o-", ms=3, label=f"{nm} ({s['p'][m].max():.0f} MPa)")
    for a, rn in zip(ax, ("fixed", "output")):
        a.set_title(f"{rn} rows, both blocks: roller loads")
        a.set_xlabel("roller (block 1, then block 2)")
        a.set_ylabel("N")
        a.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "loads.png", dpi=130)

    (HERE / "design.json").write_text(json.dumps(dict(
        wobble=asdict(w), ratio=w.ratio, gamma_deg=float(np.degrees(w.gamma)), preload_mm=w.preload,
        undercut_limit_mm=lim, flank_deg=fl, single=r1, twin=r2, reference=HPF25,
        ratio_options=ratio_table()[:40]), indent=1, default=float))
    print("wrote out/grooves.png, out/loads.png, design.json")


if __name__ == "__main__":
    main()
