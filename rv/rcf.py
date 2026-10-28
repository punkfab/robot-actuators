"""
Rolling contact fatigue of an RV reducer crankshaft's eccentric, after

  Li, Shao, Tang, Ding, Zhou, "An Investigation of the Contact Fatigue Characteristics
  of an RV Reducer Crankshaft, Considering the Hardness Gradients and Initial Residual
  Stress", Materials 15 (2022) 7850, doi:10.3390/ma15217850.

The paper's chain, rebuilt:
  hardness with depth      Thomas's two-parabola profile (their eqs 3-4): CHD, HV_surface, HV_core
  strength from hardness   yield and tensile strength from HV (eqs 5-6), Meyer index m = 2.19
  residual stress          Hertter's fit to hardness (eq 7, "ICRS"); its mirror image ("ITRS");
                           none; and a shot-peened profile read off their Fig. 11 ("SPRS")
  contact stresses         a needle roller on the eccentric as a line contact, rolled back and
                           forth over the material (reciprocating, friction 0.05, as in the paper)
  fatigue                  critical-plane search over all orientations; the modified
                           Fatemi-Socie damage (eq 13) solved for crack-initiation life N_f

One difference, stated: the paper solves the contact stresses with a 3-D elasto-plastic
FEM (ABAQUS). Here they are the closed-form elastic line-contact field with friction
(McEwen 1949, as given in Johnson, Contact Mechanics, 1985, sec. 7.1); checked against the
known surface and subsurface values in check_field(). The paper shows the
response is elastic at the rated 800 N m (Fig. 20: plastic strain ~0), so the elastic
field should be close there; it will not be at 1800 N m and above, where they find
plastic strain of 0.2-0.8%.

    python rv/rcf.py          (or: make rv-rcf)  -> out/rcf.png, out/rcf.json
"""

import json
from dataclasses import dataclass, replace
from math import pi, sqrt
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "out"
KGF = 9.80665                     # HV is kgf/mm2; x KGF -> MPa


# ---- the paper's rolling pair and material (Table 1, Section 3) -------------------------
@dataclass
class Pair:
    r_roller: float = 4.0          # mm, R1
    r_ecc: float = 16.6            # mm, R2 (eccentric cylindrical surface)
    length: float = 12.0           # mm, roller length l1
    e_roller: float = 219e3        # MPa, GCr15
    e_crank: float = 210e3         # MPa, 20CrNi2MoA
    nu_roller: float = 0.3
    nu_crank: float = 0.275

    @property
    def r_eq(self):                # 3.22 mm in the paper
        return 1 / (1 / self.r_roller + 1 / self.r_ecc)

    @property
    def e_star(self):              # plane-strain contact modulus, MPa
        return 1 / ((1 - self.nu_roller ** 2) / self.e_roller + (1 - self.nu_crank ** 2) / self.e_crank)

    def p0_of_load(self, q):
        """peak Hertz pressure (MPa) for a roller load q (N) on its full length"""
        w = q / self.length
        return sqrt(w * self.e_star / (pi * self.r_eq))

    def load_of_p0(self, p0):
        return p0 ** 2 * pi * self.r_eq / self.e_star * self.length

    def half_width(self, p0):
        """contact semi-width a, mm"""
        return 2 * self.r_eq * p0 / self.e_star


@dataclass
class Case:
    """carburized case and surface treatment of the eccentric"""
    chd: float = 1.0               # mm, depth to 550 HV
    hv_surface: float = 670.0
    hv_core: float = 450.0
    residual: str = "icrs"         # icrs | none | itrs | sprs
    m_meyer: float = 2.19


def hardness(c: Case, y):
    """HV at depth y (mm): Thomas profile, eqs 3-4 with y_HV,max = 0"""
    y = np.asarray(y, float)
    aa = (550 - c.hv_surface) / c.chd ** 2
    ba, ca = 0.0, c.hv_surface
    slope = 2 * aa * c.chd + ba
    # the second parabola has its vertex at y_cor with value HV_core (continuity of value and slope at CHD)
    d = 2 * (550 - c.hv_core) / -slope                       # y_cor - CHD
    y_cor = c.chd + d
    ab = slope / (2 * (c.chd - y_cor))
    bb = -2 * ab * y_cor
    cb = 550 - ab * c.chd ** 2 - bb * c.chd
    return np.where(y < c.chd, aa * y * y + ba * y + ca, np.where(y < y_cor, ab * y * y + bb * y + cb, c.hv_core))


def strength(c: Case, hv):
    """yield and tensile strength, MPa, from HV (eqs 5-6)"""
    s = np.asarray(hv) * KGF / 3
    k = c.m_meyer - 2
    return s * 0.1 ** k, s * (1 - k) * (12.5 * k / (1 - k)) ** k


# shot-peened residual stress, read by eye off the paper's Fig. 11 (depth mm, MPa)
SPRS_FIG11 = [(0.0, -800), (0.012, -1040), (0.025, -1170), (0.05, -1200), (0.075, -1120), (0.1, -975),
              (0.125, -825), (0.15, -700), (0.175, -595), (0.2, -500), (0.225, -430), (0.25, -375),
              (0.275, -335), (0.3, -300), (0.35, -275), (0.4, -262), (0.5, -245), (0.6, -225)]


def residual(c: Case, y):
    """residual stress along the rolling and axial directions (equal, per the paper), MPa"""
    y = np.asarray(y, float)
    dh = hardness(c, y) - c.hv_core
    icrs = np.where(dh <= 300, -1.25 * dh, 0.2857 * dh - 460)       # Hertter, eq 7
    if c.residual == "icrs":
        return icrs
    if c.residual == "itrs":
        return -icrs
    if c.residual == "none":
        return 0 * y
    if c.residual == "sprs":
        d, s = np.array(SPRS_FIG11).T
        return np.interp(y, d, s, right=float(np.interp(0.6, d, s)))
    raise ValueError(c.residual)


# ---- fatigue constants (Section 3.4) ------------------------------------------------------
B_EXP, C_EXP, K_FS = -0.087, -0.58, 1.0
CRIT = "gamma"           # critical plane: "gamma" (max shear strain amplitude, FS) or "damage" (max FD)


def fs_constants(c: Case, hv, e_mod):
    sy, sb = strength(c, hv)
    psi = np.where(sb / e_mod < 0.003, 1.0, 1.375 - 125 * sb / e_mod)
    sf, ef = 1.5 * sb, 0.59 * psi
    return dict(sy=sy, sb=sb, tau_f=sf / sqrt(3), gam_f=sqrt(3) * ef)


def life_from_damage(fd, tau_f, gam_f, g):
    """N_f from FD = tau'f/G (2N)^b + gamma'f (2N)^c (Newton in log 2N)"""
    x = np.full_like(np.asarray(fd, float), 16.0)              # x = ln(2N)
    for _ in range(80):
        f = tau_f / g * np.exp(B_EXP * x) + gam_f * np.exp(C_EXP * x) - fd
        df = tau_f / g * B_EXP * np.exp(B_EXP * x) + gam_f * C_EXP * np.exp(C_EXP * x)
        x = x - f / df
    return np.exp(x) / 2


# ---- contact stress field: Hertz line contact with sliding friction ------------------------
def line_contact(x, z, a, p0, mu):
    """stresses (MPa) at points (x along rolling, z depth >= 0, mm) under a Hertz line load of
    semi-width a and peak p0 centred at x = 0, with tangential traction mu*p along +x.
    McEwen's closed form (Johnson 1985, sec. 7.1)."""
    x, z = np.broadcast_arrays(np.asarray(x, float), np.asarray(z, float))
    A = a * a - x * x + z * z
    root = np.sqrt(A * A + 4 * x * x * z * z)
    m = np.sqrt(np.maximum(0.5 * (root + A), 0.0))
    n = np.sign(x) * np.sqrt(np.maximum(0.5 * (root - A), 0.0))
    mn = m * m + n * n
    mn = np.where(mn == 0, 1e-30, mn)
    k = p0 / a
    sx = -k * (m * (1 + (z * z + n * n) / mn) - 2 * z)
    sz = -k * m * (1 - (z * z + n * n) / mn)
    txz = k * n * (m * m - z * z) / mn
    q = -mu * k                    # McEwen's tangential terms below are for traction toward -x; flip so +mu is +x
    sx += -q * (n * (2 - (z * z - m * m) / mn) - 2 * x)
    sz += q * n * (m * m - z * z) / mn
    txz += -q * (m * (1 + (z * z + n * n) / mn) - 2 * z)
    return sx, sz, txz


def check_field():
    """known values for mu = 0: surface centre sx = sz = -p0; max principal shear 0.300 p0 at
    z = 0.786 a; max |txz| 0.250 p0 at z = 0.5 a, x = +-0.87 a"""
    a, p0 = 1.0, 1.0
    z = np.linspace(0, 2, 2001)
    sx, sz, txz = line_contact(0 * z, z, a, p0, 0.0)
    t1 = 0.5 * np.abs(sx - sz)
    xs = np.linspace(-2, 2, 801)
    X, Zg = np.meshgrid(xs, np.linspace(0.01, 1.5, 300))
    _, _, T = line_contact(X, Zg, a, p0, 0.0)
    k = np.unravel_index(np.argmax(np.abs(T)), T.shape)
    s0 = line_contact(0.0, 1e-9, a, p0, 0.0)
    return dict(surface_sx=float(s0[0]), surface_sz=float(s0[1]), tau1_max=float(t1.max()),
                tau1_depth=float(z[np.argmax(t1)]), txz_max=float(np.abs(T).max()), txz_depth=float(Zg[k]),
                txz_x=float(abs(X[k])))


# ---- critical plane --------------------------------------------------------------------------
def planes(step_deg=6.0):
    """unit normals over a hemisphere (alpha in [0, 180), beta in [0, 180]) and, for each, two
    in-plane unit vectors a, b (their eq 15)"""
    al = np.radians(np.arange(0, 180, step_deg))
    be = np.radians(np.arange(0, 180 + 1e-9, step_deg))
    A, B = np.meshgrid(al, be, indexing="ij")
    A, B = A.ravel(), B.ravel()
    n = np.stack([np.sin(B) * np.cos(A), np.sin(B) * np.sin(A), np.cos(B)], 1)
    av = np.stack([-np.sin(A), np.cos(A), 0 * A], 1)
    bv = np.stack([-np.cos(B) * np.cos(A), -np.cos(B) * np.sin(A), np.sin(B)], 1)
    return n, av, bv


def histories(pair, p0, mu, depth, travel=4.6, steps=161):
    """stress tensors over one reciprocating cycle at a material point at depth (mm): the
    roller passes forward (traction +x) then back (traction -x). Axes: x rolling, y depth,
    z axial (the paper's). Plane strain: sigma_zz = nu (sxx + syy)."""
    a = pair.half_width(p0)
    xs = np.linspace(-travel * a, travel * a, steps)
    fwd = line_contact(-xs, depth, a, p0, mu)                 # load at +xs, point at 0 -> relative x = -xs
    back = line_contact(-xs[::-1], depth, a, p0, -mu)
    sxx = np.concatenate([fwd[0], back[0]])
    syy = np.concatenate([fwd[1], back[1]])
    sxy = np.concatenate([fwd[2], back[2]])
    nu = pair.nu_crank
    szz = nu * (sxx + syy)
    S = np.zeros((len(sxx), 3, 3))
    S[:, 0, 0], S[:, 1, 1], S[:, 2, 2] = sxx, syy, szz
    S[:, 0, 1] = S[:, 1, 0] = sxy
    return S


def damage_at_depth(pair, case, p0, mu, depth, nrm, av, bv, n_theta=36):
    """modified FS damage on the most damaging plane at one depth; returns dict"""
    S = histories(pair, p0, mu, depth)
    E, nu = pair.e_crank, pair.nu_crank
    G = E / (2 * (1 + nu))
    tr = np.trace(S, axis1=1, axis2=2)
    eps = ((1 + nu) * S - nu * tr[:, None, None] * np.eye(3)) / E          # contact strains only
    rs = float(residual(case, depth))
    Sr = S.copy()
    Sr[:, 0, 0] += rs
    Sr[:, 2, 2] += rs
    # normal stress on each plane over time (with residual)
    tn = np.einsum("tij,pj->tpi", Sr, nrm)
    sig_n = np.einsum("tpi,pi->tp", tn, nrm)
    sn_max = sig_n.max(0)
    # engineering shear strain vector on each plane, projected on in-plane directions
    en = np.einsum("tij,pj->tpi", eps, nrm)
    th = np.linspace(0, pi, n_theta, endpoint=False)
    u = np.cos(th)[None, :, None] * av[:, None, :] + np.sin(th)[None, :, None] * bv[:, None, :]   # p, k, 3
    gam = 2 * np.einsum("tpi,pki->tpk", en, u)                          # gamma along u (normal part drops out: u . n = 0)
    gam_a = 0.5 * (gam.max(0) - gam.min(0)).max(1)                      # amplitude, max over in-plane direction
    fd = gam_a * (1 + K_FS * sn_max / (G * 2 * np.maximum(gam_a, 1e-12)))
    if CRIT == "damage":
        i = int(np.argmax(fd))                                          # the most damaging plane
    else:
        # Fatemi-Socie's definition: the plane of maximum shear strain amplitude (ties within
        # 0.5% broken toward the larger normal stress)
        top = gam_a >= gam_a.max() * 0.995
        i = int(np.flatnonzero(top)[np.argmax(sn_max[top])])
    return dict(depth=depth, fd=float(fd[i]), gam_a=float(gam_a[i]), sn_max=float(sn_max[i]), plane=nrm[i].tolist(),
                fd_classic=float(gam_a[i] * (1 + K_FS * sn_max[i] / float(strength(case, hardness(case, depth))[0]))),
                residual=rs, G=G)


def life(pair, case, p0, mu=0.05, depths=None, plane_step=6.0):
    """minimum crack-initiation life over depth, and the profile"""
    depths = np.arange(0.005, 0.6, 0.01) if depths is None else depths
    nrm, av, bv = planes(plane_step)
    rows = []
    for y in depths:
        d = damage_at_depth(pair, case, p0, mu, float(y), nrm, av, bv)
        hv = float(hardness(case, y))
        fc = fs_constants(case, hv, pair.e_crank)
        d["hv"] = hv
        d["nf"] = float(life_from_damage(d["fd"], float(fc["tau_f"]), float(fc["gam_f"]), d["G"]))
        rows.append(d)
    k = int(np.argmin([r["nf"] for r in rows]))
    return rows[k]["nf"], rows[k]["depth"], rows


# the paper's results (Fig. 24, Section 4.3) at rated 800 N m, P_h = 1.56 GPa
PAPER = dict(sprs=2.05e8, icrs=9.08e7, none=6.84e7, itrs=3.06e7)
PAPER_TEST_HOURS, PAPER_CRANK_RPM = 2150, 585          # 2150 h -> 7.55e7 (their equivalence)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    OUT.mkdir(exist_ok=True)
    pair = Pair()
    print(f"rolling pair: R_eq {pair.r_eq:.3f} mm (paper 3.22), E* {pair.e_star / 1e3:.1f} GPa")
    chk = check_field()
    print("field check (mu = 0, a = p0 = 1): surface sx {surface_sx:.3f} sz {surface_sz:.3f} (-1, -1); "
          "max principal shear {tau1_max:.3f} at z {tau1_depth:.3f} (0.300 at 0.786); "
          "max |txz| {txz_max:.3f} at z {txz_depth:.2f}, x {txz_x:.2f} (0.250 at 0.50, 0.87)".format(**chk))
    p0 = 1560.0
    print(f"P_h {p0 / 1e3:.2f} GPa: roller load {pair.load_of_p0(p0):.0f} N on {pair.length} mm, "
          f"half-width {pair.half_width(p0) * 1e3:.0f} um")
    c0 = Case()
    for y in (0.0, 0.5, 1.0, 1.5, 2.0):
        sy, sb = strength(c0, hardness(c0, y))
        print(f"  depth {y:.1f} mm: HV {float(hardness(c0, y)):5.0f}, yield {float(sy):5.0f} MPa, tensile {float(sb):5.0f} MPa, "
              f"residual (Hertter) {float(residual(c0, y)):6.0f} MPa")

    print(f"\nlife at P_h = {p0 / 1e3:.2f} GPa, friction 0.05 (paper Fig. 24 alongside):")
    res = {}
    for rs in ("sprs", "icrs", "none", "itrs"):
        nf, y, rows = life(pair, replace(c0, residual=rs), p0)
        res[rs] = dict(nf=nf, depth=y, rows=rows)
        print(f"  {rs:5}: N_f {nf:.3g} at {y * 1e3:.0f} um   paper {PAPER[rs]:.3g}   ({nf / PAPER[rs] - 1:+.0%})")
    test = PAPER_TEST_HOURS * 60 * PAPER_CRANK_RPM
    print(f"  the paper's test: {PAPER_TEST_HOURS} h at {PAPER_CRANK_RPM} rpm = {test:.3g} crank revolutions")

    print("\nfriction, ICRS:")
    fr = {}
    for mu in (0.0, 0.05, 0.2):
        nf, y, _ = life(pair, c0, p0, mu=mu)
        fr[mu] = (nf, y)
        print(f"  mu {mu:4.2f}: N_f {nf:.3g} at {y * 1e3:.0f} um")

    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    yy = np.linspace(0, 2.2, 300)
    ax[0].plot(yy, hardness(c0, yy), color="#1f4e79")
    ax[0].axhline(550, color="0.6", lw=0.7, ls="--")
    ax[0].set(xlabel="depth (mm)", ylabel="HV", title="hardness (Thomas, CHD 1.0 mm)")
    y6 = np.linspace(0, 0.6, 200)
    for rs, col in (("sprs", "#b03a2e"), ("icrs", "#c87a2f"), ("none", "0.5"), ("itrs", "#1f4e79")):
        ax[1].plot(y6, residual(replace(c0, residual=rs), y6), color=col, label=rs.upper())
        rows = res[rs]["rows"]
        ax[2].semilogy([r["depth"] * 1e3 for r in rows], [r["nf"] for r in rows], color=col, label=rs.upper())
        ax[2].plot([res[rs]["depth"] * 1e3], [PAPER[rs]], "o", color=col, mfc="none")
    ax[1].set(xlabel="depth (mm)", ylabel="residual stress (MPa)", title="initial residual stress")
    ax[1].legend(frameon=False, fontsize=8)
    ax[2].set(xlabel="depth (um)", ylabel="crack-initiation life N_f", title="life by depth (o = paper's minimum)")
    ax[2].legend(frameon=False, fontsize=8)
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUT / "rcf.png", dpi=150)
    (OUT / "rcf.json").write_text(json.dumps(dict(
        check=chk, p0=p0, lives={k: dict(nf=v["nf"], depth=v["depth"], paper=PAPER[k]) for k, v in res.items()},
        friction={str(k): v for k, v in fr.items()}), indent=1))
    print(f"\nwrote {OUT / 'rcf.png'}, {OUT / 'rcf.json'}")


if __name__ == "__main__":
    main()
