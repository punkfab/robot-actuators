"""
wobble2: a single wobble ring between two rows of off-the-shelf rollers, a nutating reducer.

    input shaft -> carrier with 2 ball bearings (radial axles) rolling on the ring's FRONT wave
    wobble ring  -> a band cut from a sphere, seated in a spherical bore, kept from spinning by two
                    pins in axial slots; its BACK end face is a second wave
    output       -> carrier with Z_o ball bearings rolling on the back wave

Kinematics. A rigid, non-spinning ring tilted theta can only present a first-harmonic height
change, R theta cos(alpha - psi), to anything touching it. A wave of N lobes, -R theta cos(N alpha),
cancels that for a roller at alpha when psi = (1 - N) alpha (or (1 + N) alpha), and the same
holds for every roller spaced 2 pi / (N -+ 1) apart. So Z rollers on a wave of N = Z +- 1 lobes
nutate the ring Z times per carrier turn:

    psi = s_i Z_i phi = s_o Z_o phi_o        ratio  i = phi / phi_o = (s_o Z_o) / (s_i Z_i)

with s = -1 for N = Z + 1 and +1 for N = Z - 1. The ratio is a roller count, like a gear pair
(Z_o / Z_i), and needs the two waves to have DIFFERENT lobe counts; the clocking between the
faces (back_clock, 30 deg by default, after the sketch) doesn't change it. The input row needs
three bearings (a 4-lobe front wave): two are either collinear (dead centres) or lopsided.

The ring must not spin: its torque balance would otherwise force input torque = output torque.
The pins carry the reaction, and the spherical seat locates it; the rollers squeeze it from both
sides, which preloads both rows and takes out backlash. Rollers are radial (axes through the
nutation centre O), so the faces are ruled through O and the analysis is done on the rollers'
pitch sphere, radius R (the bearing's mid-plane). Bearings are cylinders, not cones: they slip
by about +-B/(2R) across their width, which this model ignores.

    python wobble2/design.py          (or: make wobble2)   -> out/*.png, design.json
    python wobble2/design.py --sweep  also re-run the options sweep (~3 min)
"""

import json
import sys
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"

E_STAR = 210e3 / (2 * (1 - 0.3 ** 2))          # MPa, steel on steel
P_DYN = 1500.0                                  # MPa, long-life contact pressure on the hardened faces


@dataclass
class Bearing:
    name: str
    d: float
    D: float
    B: float
    C: float          # N, basic dynamic load rating
    C0: float         # N, basic static load rating
    src: str = ""


# SKF ratings as listed by distributors (check the exact suffix before relying on them)
B624 = Bearing("624-2Z", 4, 13, 5, 936, 290, "SKF 624-Z via im-tek / simplybearings listings")
B626 = Bearing("626-2Z", 6, 19, 6, 2340, 950, "SKF 626-2Z/C3 via RS listing")
B623 = Bearing("623-2Z", 3, 10, 4, 540, 180, "SKF 623-2Z via coolenbearings / qualitybearingsonline listings")
B693 = Bearing("693-2Z", 3, 8, 4, 559, 179, "FAG 693-2Z via sp-spareparts listing")
# stud-type cam followers (needle-roller track rollers on a stud), IKO; outer ring width per the CF series table
CF4R = Bearing("IKO CF4R cam follower (crowned)", 4, 12, 8, 2070, 1590, "IKO CF4R via acorn-ind listing")
CF3 = Bearing("IKO CF3 cam follower", 3, 10, 7, 1500, 1020, "IKO CF3B via qualitybearingsonline listing")
B6000 = Bearing("6000-2RS", 10, 26, 8, 4750, 1960, "SKF 6000-2RSL via RS listing")


@dataclass
class Row:
    slots: int          # Z: roller positions the wave is conjugate to, 2 pi / Z apart
    lobes: int          # N = Z +- 1
    elev: float         # deg, roller axes on this cone (input below the equator, output above)
    bearing: Bearing
    side: int           # +1: the face is below the roller (output), -1: above it (input)
    present: tuple = () # which slots carry a bearing (empty: all)

    @property
    def s(self):
        if self.lobes == self.slots + 1:
            return -1
        if self.lobes == self.slots - 1:
            return 1
        raise ValueError("lobes must be slots +- 1")

    @property
    def idx(self):
        return list(self.present) if self.present else list(range(self.slots))

    @property
    def rollers(self):
        return len(self.idx)


@dataclass
class Wobble2:
    R: float = 55.0                     # mm pitch sphere (bearing mid-planes)
    theta: float = 1.8                  # deg nutation angle
    # three input bearings on a 4-lobe front wave. Two don't work: 180 deg apart (3-lobe wave) they
    # can't resist the ring tilting about the line through them, and every 60 deg of input that tilt
    # lines up with the output's rolling (a dead centre); 120 deg apart they push the ring on one
    # side only and half the output rollers unload. See the conditioning printout.
    inp: Row = field(default_factory=lambda: Row(3, 4, -24.0, B6000, -1))
    # output: stud-type cam followers (needle rollers, crowned), the bearing-on-a-stud of the sketch
    out: Row = field(default_factory=lambda: Row(12, 11, 24.0, CF4R, 1))
    back_clock: float = 30.0            # deg, the output carrier's angle at phi = 0 (clocks the back wave)
    preload_N: float = 60.0             # N per output roller from the squeeze, no load
    pin_r: float = 0.0                  # mm, anti-rotation pin radius from the axis (0: = seat radius)
    mu_roll: float = 0.0015             # ball-bearing rollers, rolling friction per unit load
    mu_pin: float = 0.12                # pin sliding in its slot (steel on steel, greased)

    @property
    def ratio(self):
        return (self.out.s * self.out.slots) / (self.inp.s * self.inp.slots)

    def gamma(self, row):
        return np.arcsin(row.bearing.D / 2 / self.R)


# ---------------------------------------------------------------- kinematics

def Rz(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])


def Rx(a):
    c, s = np.cos(a), np.sin(a)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def psi(w, phi):
    return w.inp.s * w.inp.slots * phi


def ring_rot(w, phi):
    """the wobble ring: tilted theta, tilt direction psi, no spin"""
    p = psi(w, phi)
    return Rz(p) @ Rx(np.radians(w.theta)) @ Rz(-p)


def out_angle(w, phi):
    return phi / w.ratio + np.radians(w.back_clock)


def carrier_rot(w, row, phi):
    return Rz(phi) if row is w.inp else Rz(out_angle(w, phi))


def axes(row):
    b = 2 * np.pi * np.array(row.idx) / row.slots
    e = np.radians(row.elev)
    return np.stack([np.cos(b) * np.cos(e), np.sin(b) * np.cos(e), np.full_like(b, np.sin(e))], 1)


def path(w, row, phis, u):
    """roller axis in the RING frame"""
    return np.array([ring_rot(w, p).T @ carrier_rot(w, row, p) @ u for p in np.atleast_1d(phis)])


def frame(w, row, phi, U, h=1e-5):
    """for roller axes U (n,3) at phi, in the ring frame: p, unit tangent t, n_up (sphere-tangent
    normal toward higher elevation), curvature k toward n_up (1/rad), speed |dp/dphi| (rad/rad)"""
    P = np.stack([(ring_rot(w, q).T @ carrier_rot(w, row, q)) @ U.T for q in (phi - h, phi, phi + h)])
    p = P[1]
    d1 = (P[2] - P[0]) / (2 * h)
    d2 = (P[2] - 2 * P[1] + P[0]) / h ** 2
    sp = np.sqrt((d1 * d1).sum(0))
    t = d1 / sp
    n = np.cross(p.T, t.T).T
    up = np.array([[0], [0], [1.0]]) - p * p[2]
    n = n * np.where((n * up).sum(0) < 0, -1.0, 1.0)
    k = (d2 * n).sum(0) / sp ** 2
    return p.T, t.T, n.T, k, sp


def face(w, row, samples=6000):
    """the ring face this row rolls on, ring frame: (azimuth, elevation) in rad along roller 0's
    path, plus the tightest convex radius (mm, at R) a roller meets"""
    u = axes(row)[:1]
    span = 2 * np.pi * (abs(w.ratio) if row is w.out else 1)          # roller 0 goes once round the ring
    phis = np.linspace(0, span, samples, endpoint=False)
    g = w.gamma(row)
    pts, kk = [], []
    for q in phis:
        p, t, n, k, sp = frame(w, row, q, u)
        pts.append(np.cos(g) * p[0] - row.side * np.sin(g) * n[0])
        kk.append(k[0] * row.side)
    pts = np.array(pts)
    az = np.unwrap(np.arctan2(pts[:, 1], pts[:, 0]))
    el = np.arcsin(np.clip(pts[:, 2], -1, 1))
    return az, el, np.array(kk)


def check_face(w, row):
    """every roller of the row must trace the same face (the conjugacy check), and the face must
    not cusp: the path curvature radius has to exceed the roller radius"""
    U = axes(row)
    n = 4000
    span = 2 * np.pi * (abs(w.ratio) if row is w.out else 1)
    phis = np.linspace(0, span, n, endpoint=False)
    P0 = np.array([frame(w, row, q, U[:1])[0][0] for q in phis])
    from scipy.optimize import minimize_scalar
    dphi = phis[1] - phis[0]
    worst = 0.0
    for j in range(1, len(U)):
        for q in phis[::200]:
            pj = frame(w, row, q, U[j:j + 1])[0][0]
            i0 = int(np.argmin(np.linalg.norm(P0 - pj, axis=1)))
            r = minimize_scalar(lambda a: np.linalg.norm(frame(w, row, a, U[:1])[0][0] - pj),
                                bounds=(phis[i0] - 2 * dphi, phis[i0] + 2 * dphi), method="bounded",
                                options=dict(xatol=1e-12))
            worst = max(worst, r.fun)
    ks = np.abs([frame(w, row, q, U[:1])[3][0] for q in phis[::4]])
    lim = w.R / ks.max()
    slope = []
    for q in phis[::8]:
        p, t, nn, k, sp = frame(w, row, q, U[:1])
        az = np.array([-p[0, 1], p[0, 0], 0]) / np.hypot(p[0, 0], p[0, 1])
        slope.append(np.degrees(np.arctan2(abs(t[0, 2]), abs(t[0] @ az))))
    return dict(conjugacy_um=worst * w.R * 1e3, undercut_r_mm=lim, roller_r_mm=row.bearing.D / 2,
                steepest_deg=max(slope))


# ---------------------------------------------------------------- contacts and statics

def k_contact(row):
    """ball-bearing roller on a face: Palmgren line contact over the bearing width, plus the
    bearing's own radial compliance (approx. as a second contact in series)"""
    return (row.bearing.B ** 0.8 / (2 * 3.84e-5)) ** (1 / 0.9)


def contacts(w, phi):
    """all roller contacts, world frame (ring at its ideal pose for phi). Arrays: C (n,3) mm,
    n (n,3) unit normal the face pushes the roller along, t path tangent, k curvature (+ conforming),
    v roller travel along the face per input radian (mm), row index (0 input, 1 output)"""
    out = {k: [] for k in ("C", "n", "t", "k", "v", "row")}
    Rr = ring_rot(w, phi)
    for ri, row in enumerate((w.inp, w.out)):
        p, t, n, k, sp = frame(w, row, phi, axes(row))
        g = w.gamma(row)
        s = row.side
        nin = s * n                        # face below the roller (s=+1) pushes it up, above (s=-1) down
        C = w.R * (np.cos(g) * p - s * np.sin(g) * n)
        out["C"].append((Rr @ C.T).T)
        out["n"].append((Rr @ nin.T).T)
        out["t"].append((Rr @ t.T).T)
        out["k"].append(k * s)
        out["v"].append(w.R * sp)
        out["row"].append(np.full(len(k), ri))
    return {k: np.concatenate(v) for k, v in out.items()}


def hertz(w, Q, k, rows):
    """peak line-contact pressure, MPa: bearing outer race (radius D/2) on the face"""
    p = np.zeros_like(Q)
    for ri, row in enumerate((w.inp, w.out)):
        m = rows == ri
        r = row.bearing.D / 2
        g = w.gamma(row)
        kk = k[m]
        wall = (1 / np.maximum(np.abs(kk), 1e-9) + np.sign(kk) * g) * w.R * np.sign(kk)
        inv = np.maximum(1 / r - 1 / wall, 1e-6)
        p[m] = np.sqrt(np.maximum(Q[m], 0) * E_STAR * inv / (np.pi * row.bearing.B))
    return p


_PRE = {}


def preload_mm(w):
    """the squeeze interference (same at every roller) that leaves the output rollers averaging
    preload_N once the ring has settled between the two rows (found by bisection, cached)"""
    key = (w.R, w.theta, w.inp.slots, w.inp.lobes, tuple(w.inp.idx), w.inp.elev, w.inp.bearing.name, w.out.slots,
           w.out.lobes, w.out.elev, w.out.bearing.name, w.preload_N, w.back_clock)
    if key in _PRE:
        return _PRE[key]
    lo, hi = 0.0, 0.05
    for _ in range(40):
        mid = (lo + hi) / 2
        Q = solve(w, 0.0, pre=mid)["Q"]
        mean = Q[-w.out.rollers:].mean()
        lo, hi = (mid, hi) if mean < w.preload_N else (lo, mid)
    _PRE[key] = (lo + hi) / 2
    return _PRE[key]


def solve(w, phi=0.0, T=0.0, Fa=0.0, k_seat=1e6, k_pin=1e9, k_ob=1e7, k_ax=50.0, pre=None):
    """statics with the input locked. Unknowns: ring (v, r) and output carrier (v, r), N mm rad.
    The ring sits in a spherical seat (translations held, k_seat N/mm) and two pins hold its spin
    (k_pin N mm/rad about z); the output carrier is on its main bearing (k_ob on everything but its
    rotation). Load: output torque T (N mm) and axial Fa (N) on the output carrier.
    Every roller starts at the same interference, set so the output rollers carry preload_N."""
    c = contacts(w, phi)
    K = np.array([k_contact(w.inp) if r == 0 else k_contact(w.out) for r in c["row"]])
    pre = preload_mm(w) if pre is None else pre
    n, C = c["n"], c["C"]
    # delta = pre - (d_roller - d_ring).n_in ; roller on input (fixed) or output (dofs 6..11), face on ring (0..5)
    J = np.zeros((len(K), 12))
    rc = np.cross(C, n)
    # d_rel = d_roller - d_ring, delta = pre - d_rel . n_in: +[n, C x n] on the ring's dofs,
    # -[n, C x n] on the output carrier's (the input carrier is held)
    J[:, 0:3], J[:, 3:6] = n, rc
    out = c["row"] == 1
    J[out, 6:9], J[out, 9:12] = -n[out], -rc[out]
    # anti-rotation: a pure nutation (tilt theta, tilt direction turning, no spin) never rotates about
    # the bisector b of the housing axis and the ring axis (the constant-velocity-joint plane), so
    # that is the one rotation the pins lock. (Two pins in straight axial slots lock world z instead,
    # a Cardan joint, which fights the nutation at second order in theta: ~17 um at the faces here.)
    a = ring_rot(w, phi) @ np.array([0, 0, 1.0])
    b = (a + np.array([0, 0, 1.0]))
    b /= np.linalg.norm(b)
    # the seat locates the ring radially; it has axial clearance, so the two rows squeeze the ring
    # against each other and share the preload
    Kc = np.diag([k_seat] * 2 + [k_ax] + [0, 0, 0] + [k_ob] * 5 + [0])
    Kc[3:6, 3:6] = k_pin * np.outer(b, b)
    G = np.zeros(12)
    G[8], G[11] = Fa, T
    def energy(x):
        dp = np.maximum(pre + J @ x, 0)
        return (K * 9 / 19 * dp ** (19 / 9)).sum() + 0.5 * x @ Kc @ x - G @ x

    x = np.zeros(12)
    U = energy(x)
    for _ in range(300):
        dp = np.maximum(pre + J @ x, 0)
        Q = K * dp ** (10 / 9)
        g = J.T @ Q + Kc @ x - G
        H = (J * (K * 10 / 9 * dp ** (1 / 9))[:, None]).T @ J + Kc + 1e-6 * np.eye(12)
        dx = -np.linalg.solve(H, g)
        step = 1.0
        while True:                                            # backtrack on the energy
            U2 = energy(x + step * dx)
            if U2 <= U + 1e-4 * step * (g @ dx) or step < 1e-10:
                break
            step /= 2
        x, U = x + step * dx, U2
        if np.abs(step * dx).max() < 1e-13:
            break
    g = J.T @ (K * np.maximum(pre + J @ x, 0) ** (10 / 9)) + Kc @ x - G
    d = pre + J @ x
    Q = K * np.maximum(d, 0) ** (10 / 9)
    F = Q[:, None] * n                                            # on each roller
    inp = c["row"] == 0
    tau_in = -np.cross(C[inp], F[inp]).sum(0)[2]                  # what the input must supply about z
    # the pins' torque on the ring (along the bisector)
    M_pin = -k_pin * (b @ x[3:6])
    p = hertz(w, Q, c["k"], c["row"])
    return dict(x=x, Q=Q, p=p, row=c["row"], tau_in=tau_in, M_pin=M_pin, c=c, ok=bool(np.abs(g).max() < 1e-3 * (1 + abs(T))))


def worst(w, phases=16, **load):
    out = None
    for ph in np.linspace(0, 2 * np.pi / w.inp.slots, phases, endpoint=False):
        s = solve(w, ph, **load)
        if out is None or s["Q"].max() > out["Q"].max():
            out = s
    return out


def bearing_life_h(b, P, rpm):
    if P <= 0:
        return np.inf
    return (b.C / P) ** 3 * 1e6 / (60 * rpm)


def rating(w, T, n_in=3000.0, phases=16):
    """at output torque T (N m) and input speed n_in (rpm): peak loads, pressure, L10 of each row's
    bearings (cubic mean load over a cycle), pin force"""
    Qs, ps = [], []
    for ph in np.linspace(0, 2 * np.pi / w.inp.slots, phases, endpoint=False):
        s = solve(w, ph, T=T * 1e3)
        Qs.append(s["Q"])
        ps.append(s["p"])
    Qs, ps = np.array(Qs), np.array(ps)
    rows = s["row"]
    res = {}
    for ri, row in enumerate((w.inp, w.out)):
        m = rows == ri
        Q = Qs[:, m]
        Pm = np.cbrt((Q ** 3).mean())                             # every roller passes every phase
        rpm_roller = (n_in if ri == 0 else n_in / abs(w.ratio)) * w.R / (row.bearing.D / 2)
        res[["input", "output"][ri]] = dict(Qmax=float(Q.max()), Pmean=float(Pm), p_max=float(ps[:, m].max()),
                                            C0_ratio=float(Q.max() / row.bearing.C0), rpm=rpm_roller,
                                            L10_h=float(bearing_life_h(row.bearing, Pm, rpm_roller)))
    return res


def capacity(w, n_in=3000.0):
    """output torque at which the first limit is reached: a roller at its bearing's C0, or the
    face pressure at P_DYN"""
    def lim(T):
        r = rating(w, T, n_in, phases=8)
        return max(max(v["C0_ratio"] for v in r.values()), max(v["p_max"] for v in r.values()) / P_DYN)
    lo, hi = 0.0, 1.0
    while lim(hi) < 1:
        lo, hi = hi, hi * 2
        if hi > 1e4:
            return np.inf, "none"
    for _ in range(16):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if lim(mid) < 1 else (lo, mid)
    r = rating(w, lo, n_in, phases=8)
    who = max(((k, max(v["C0_ratio"], v["p_max"] / P_DYN)) for k, v in r.items()), key=lambda kv: kv[1])[0]
    gov = "C0" if r[who]["C0_ratio"] >= r[who]["p_max"] / P_DYN else "face pressure"
    return lo, f"{who} rollers, {gov}"


def efficiency(w, T, phases=16):
    """forward efficiency at output torque T (N m): rolling loss of every roller plus the pins
    sliding in their slots (they carry the ring's reaction torque and ride its rocking)"""
    th = np.radians(w.theta)
    Rp = w.pin_r or seat_radius(w)
    effs = []
    for ph in np.linspace(0, 2 * np.pi / w.inp.slots, phases, endpoint=False):
        s = solve(w, ph, T=T * 1e3)
        loss_roll = w.mu_roll * (s["Q"] * s["c"]["v"]).sum()
        # each pin sits at the ring's equator; its axial speed per input radian, averaged over a turn,
        # is (2/pi) Rp theta |dpsi/dphi|; two pins share the reaction torque
        F_pin = abs(s["M_pin"]) / (2 * Rp)
        loss_pin = w.mu_pin * 2 * F_pin * (2 / np.pi) * Rp * th * abs(w.inp.s * w.inp.slots)
        work = T * 1e3 / abs(w.ratio)
        effs.append(work / (work + loss_roll + loss_pin))
    return float(np.mean(effs))


def conditioning(w, phi):
    """how well the rollers hold the drive at phi: the smallest singular value of the contact
    Jacobian over the ring's rotations and the output's rotation (pins on the bisector added).
    Zero means a mechanism with the input locked: a dead centre."""
    c = contacts(w, phi)
    n, C = c["n"], c["C"]
    rc = np.cross(C, n)
    J = np.zeros((len(n), 4))
    J[:, 0:3] = rc
    o = c["row"] == 1
    J[o, 3] = -rc[o, 2]
    a = ring_rot(w, phi) @ np.array([0, 0, 1.0])
    b = (a + np.array([0, 0, 1.0])) / np.linalg.norm(a + np.array([0, 0, 1.0]))
    s = np.linalg.svd(np.vstack([J, np.r_[b * 100, 0]]), compute_uv=False)
    return s[-1]


def seat_radius(w):
    """the ring's outer sphere (and the pins' radius), just beyond the bearings' outer ends"""
    return w.R + max(w.inp.bearing.B, w.out.bearing.B) / 2 + 4.0


# ---------------------------------------------------------------- main

def sweep():
    rows = []
    for R in (45.0, 55.0):
        for th in (1.5, 2.0, 2.5, 3.0):
            for Zo in (8, 10, 12, 14, 16, 18, 20, 24):
                for No in (Zo - 1, Zo + 1):
                    for brg in (B693, CF3, CF4R):
                        if Zo * (brg.D + 3) > 2 * np.pi * R:          # bearings must fit round the carrier
                            continue
                        w = Wobble2(R=R, theta=th, out=Row(Zo, No, 24.0, brg, 1))
                        co = check_face(w, w.out)
                        if co["undercut_r_mm"] < 1.25 * brg.D / 2:
                            continue
                        if co["steepest_deg"] < 18:
                            continue
                        rows.append((R, th, Zo, No, brg.name, w.ratio, co["steepest_deg"], co["undercut_r_mm"]))
    return rows


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(exist_ok=True)
    w = Wobble2()
    print(f"ratio {w.ratio:+.2f}  (input {w.inp.rollers} rollers in {w.inp.slots} slots on a {w.inp.lobes}-lobe wave, output "
          f"{w.out.rollers} on {w.out.lobes}); ring nutates {abs(w.inp.s * w.inp.slots)}x per input turn; theta {w.theta} deg, R {w.R} mm")
    cond = {}
    for nm, ww in (("2 bearings 180 deg apart (3-lobe)", replace(w, inp=Row(2, 3, w.inp.elev, w.inp.bearing, -1))),
                   ("2 bearings 120 deg apart (4-lobe)", replace(w, inp=Row(3, 4, w.inp.elev, w.inp.bearing, -1, present=(0, 1)))),
                   ("3 bearings 120 deg apart (4-lobe)", w)):
        v = [conditioning(ww, p) for p in np.linspace(0, 2 * np.pi, 361)]
        Qo = solve(ww, 0.4)["Q"][-ww.out.rollers:]
        cond[nm] = dict(conditioning=float(min(v) / max(v)), output_unloaded_at_preload=int((Qo < 0.05 * Qo.max()).sum()))
        print(f"  {nm:34s}: drive conditioning worst/best {cond[nm]['conditioning']:.3f}"
              + ("  (dead centres)" if cond[nm]["conditioning"] < 0.01 else "")
              + f"; output rollers unloaded at preload {cond[nm]['output_unloaded_at_preload']}/{ww.out.rollers}")
    chk = {}
    for nm, row in (("input", w.inp), ("output", w.out)):
        chk[nm] = check_face(w, row)
        c = chk[nm]
        print(f"  {nm:6s} face: rollers on one face within {c['conjugacy_um']:.1f} um; steepest flank {c['steepest_deg']:.1f} deg;"
              f" roller r {c['roller_r_mm']} mm vs undercut limit {c['undercut_r_mm']:.1f} mm")
    s = solve(w, 0.3, T=20e3)
    print(f"  torque check: 20 N m out needs {abs(s['tau_in']) / 1e3:.4f} N m in, expected {20 / abs(w.ratio):.4f}; "
          f"pins carry {abs(s['M_pin']) / 1e3:.2f} N m")
    s0 = worst(w)
    print(f"  preload only: input rollers {s0['Q'][s0['row'] == 0].max():.0f} N, output {s0['Q'][s0['row'] == 1].min():.0f}-"
          f"{s0['Q'][s0['row'] == 1].max():.0f} N")
    cap, gov = capacity(w)
    print(f"  torque capacity {cap:.1f} N m (limit: {gov})")
    rat = {}
    for T in (10.0, cap * 0.5):
        r = rating(w, T)
        rat[f"{T:.1f}"] = r
        print(f"  at {T:5.1f} N m, 3000 rpm in: " + "; ".join(
            f"{k} rollers Qmax {v['Qmax']:.0f} N, p {v['p_max']:.0f} MPa, {v['rpm']:.0f} rpm, L10 {v['L10_h']:.0f} h"
            for k, v in r.items()))
    eff = {f"{T} N m": efficiency(w, T) for T in (2, 5, 10, 20, 30)}
    print("  efficiency: " + ", ".join(f"{k} {v * 100:.1f}%" for k, v in eff.items()))
    sw = sweep() if "--sweep" in sys.argv else json.loads((HERE / "design.json").read_text()).get("options", []) \
        if (HERE / "design.json").exists() else []
    print("  options (no cusp at 1.25x the roller radius): R, theta, Z_o, lobes, bearing -> ratio, flank")
    for r in sw:
        print(f"    R {r[0]:.0f}  theta {r[1]}  Z_o {r[2]:2d}/{r[3]:2d}  {r[4]}  -> {r[5]:+5.1f}:1, flank {r[6]:.0f} deg")

    # figure: both faces unrolled on the pitch sphere, with the rollers at phi = 0
    fig, ax = plt.subplots(figsize=(11, 4.2))
    for row, col, nm in ((w.inp, "#1c7ed6", "front face (input rollers)"), (w.out, "#d9480f", "back face (output rollers)")):
        az, el, _ = face(w, row)
        a = np.mod(np.degrees(az), 360)
        o = np.argsort(a)
        ax.plot(a[o], np.degrees(el[o]), ".", ms=1, color=col, label=nm)
        p, *_ = frame(w, row, 0.0, axes(row))
        ax.scatter(np.mod(np.degrees(np.arctan2(p[:, 1], p[:, 0])), 360), np.degrees(np.arcsin(p[:, 2])), s=40,
                   color=col, edgecolors="k", zorder=3)
    ax.set_xlabel("azimuth on the ring (deg)")
    ax.set_ylabel("elevation (deg)")
    ax.set_title(f"wobble2 ring faces: {w.inp.lobes}-lobe front, {w.out.lobes}-lobe back (ring frame, phi = 0)")
    ax.legend(loc="center right", fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "faces.png", dpi=130)

    (HERE / "design.json").write_text(json.dumps(dict(
        wobble2=asdict(w), ratio=w.ratio, faces=chk, conditioning=cond, capacity_Nm=cap, governs=gov, ratings=rat, efficiency=eff,
        seat_radius=seat_radius(w), options=sw), indent=1, default=float))
    print("wrote out/faces.png, design.json")


if __name__ == "__main__":
    main()
