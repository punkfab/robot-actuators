"""
Nonlinear 2-D field solution for the liquid-cooled QDD motor: where does the iron
actually saturate, and what does trading copper for iron buy?

design.py assumed a saturation knee (an air-gap shear of 90 kPa). This replaces the
assumption: it meshes the cross-section of a design.Design, solves the magnetostatic
field with a real B-H curve, and returns torque against current.

  geometry    from design.Design: both topologies (stator outside / outrunner), tooth
              tips, 24s/22p double-layer concentrated winding (one coil per tooth)
  mesh        gmsh, fragment + weld (the recipe from motor/fea.py)
  solve       A_z, P1 triangles, Newton on nu(B^2), assembled directly in numpy
  iron        Brauer fit nu = k1 exp(k2 B^2) + k3, continued at the slope of air past
              saturation [APPROX: typical non-oriented electrical steel]
  torque      Arkkio's method over the air-gap elements

Limits: 2-D (no end leakage: a 25 mm stack on a 71 mm bore will lose some torque to
it), one rotor position (no ripple average), magnets at a fixed temperature, no
demagnetization check, no iron loss.

    python qdd-liquid/fea.py            (or: make qdd-liquid-fea; a few minutes)
"""

import json
import sys
from dataclasses import replace
from math import cos, pi, sin, sqrt, tau
from pathlib import Path

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import spsolve

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402

MU0 = 4e-7 * pi
NU0 = 1 / MU0
AIR, IRON_S, IRON_R, COIL, GAP, MAGNET = range(6)

# Brauer coefficients for H = nu(B^2) B. [APPROX] 1.5 T at ~1.3 kA/m, 1.8 T at ~8 kA/m.
STEELS = {
    "silicon steel": dict(k1=3.8, k2=2.17, k3=396.2, scale=1.0),
    # cobalt-iron saturates ~15% higher; modelled as the same curve with B stretched
    "cobalt iron": dict(k1=3.8, k2=2.17, k3=396.2, scale=1.15),
}


def _bh(steel):
    """Tables of nu and d nu / d(B^2) against B^2, with the curve continued at
    dH/dB = nu0 beyond the point where the Brauer fit gets that steep."""
    k1, k2, k3, sc = (STEELS[steel][k] for k in ("k1", "k2", "k3", "scale"))
    b = np.linspace(1e-4, 6.0, 6000)
    bs = b / sc
    h = (k1 * np.exp(np.minimum(k2 * bs * bs, 50)) + k3) * bs
    dh = np.gradient(h, b)
    i = int(np.argmax(dh >= NU0))
    h[i:] = h[i] + NU0 * (b[i:] - b[i])
    nu = h / b
    return b * b, nu, np.gradient(nu, b * b), b[i]


class Motor:
    """Mesh and material layout for one design."""

    def __init__(self, d: dz.Design, steel="silicon steel", br=1.17, mur=1.05,
                 open_frac=0.2, magnet_arc=0.85, h_gap=0.25, h_coarse=1.6):
        self.d, self.steel, self.br, self.mur = d, steel, br, mur
        self.open_frac, self.magnet_arc = open_frac, magnet_arc
        rg = d.r_gap * 1e3
        s_in, s_out = (x * 1e3 for x in d.slot_radii)
        if d.outrunner:                  # hub, yoke, slots, tips | gap | magnets, rotor yoke
            self.r = dict(stator_a=s_in - d.stator_yoke, slot_a=s_in, slot_b=s_out,
                          bore=s_out + d.tooth_tip, mag_a=rg, mag_b=rg + d.magnet,
                          outer=rg + d.magnet + d.rotor_yoke)
            self.gap = (self.r["bore"], rg)
        else:                            # rotor yoke, magnets | gap | tips, slots, yoke
            self.r = dict(inner=rg - d.magnet - d.rotor_yoke, mag_a=rg - d.magnet, mag_b=rg,
                          bore=s_in - d.tooth_tip, slot_a=s_in, slot_b=s_out,
                          outer=s_out + d.stator_yoke)
            self.gap = (rg, self.r["bore"])
        assert not (d.slot_w and d.outrunner), "parallel slots: stator-outside only"
        self.open_w = 2 * self.r["bore"] * sin(open_frac * pi / d.slots)      # slot opening, mm
        self._mesh(h_gap, h_coarse)
        self._materials()
        self._winding()

    # ---- regions -------------------------------------------------------------
    def classify(self, x, y):
        d, R = self.d, self.r
        r, th = np.hypot(x, y), np.arctan2(y, x)
        sp_, pp_ = tau / d.slots, tau / d.poles
        in_arc = lambda pitch, frac: np.abs((th + pitch / 2) % pitch - pitch / 2) <= 0.5 * frac * pitch
        if d.slot_w:                     # parallel-sided: local (along, across) the slot axis
            dth = (th + sp_ / 2) % sp_ - sp_ / 2
            u, v = r * np.cos(dth), r * np.sin(dth)
            slot = (u > R["slot_a"]) & (u < R["slot_b"]) & (np.abs(v) <= d.slot_w / 2)
        else:
            slot = (r > R["slot_a"]) & (r < R["slot_b"]) & in_arc(sp_, d.slot_frac)
        tip_lo, tip_hi = sorted((R["bore"], R["slot_b"] if d.outrunner else R["slot_a"]))
        if d.slot_w:
            opening = (r > R["bore"]) & (u <= R["slot_a"]) & (np.abs(v) <= self.open_w / 2)
        else:
            opening = (r > tip_lo) & (r < tip_hi) & in_arc(sp_, self.open_frac)
        magnet = (r > R["mag_a"]) & (r < R["mag_b"]) & in_arc(pp_, self.magnet_arc)
        reg = np.full(x.shape, AIR)
        if d.outrunner:
            stator = (r > R["stator_a"]) & (r < R["bore"])
            rotor = (r > R["mag_b"])
        else:
            stator = (r > R["bore"])
            rotor = (r > R["inner"]) & (r < R["mag_a"])
        reg[stator] = IRON_S
        reg[rotor] = IRON_R
        reg[opening] = AIR
        reg[slot] = COIL
        reg[(r > self.gap[0]) & (r < self.gap[1])] = GAP
        reg[magnet] = MAGNET
        return reg

    def _mesh(self, h_gap, h_coarse):
        import gmsh
        d, R, mm = self.d, self.r, 1e-3
        gmsh.initialize()
        try:
            gmsh.option.setNumber("General.Terminal", 0)
            occ = gmsh.model.occ

            def annulus(r_in, r_out):
                o = occ.addDisk(0, 0, 0, r_out * mm, r_out * mm)
                if r_in <= 0:
                    return o
                i = occ.addDisk(0, 0, 0, r_in * mm, r_in * mm)
                return occ.cut([(2, o)], [(2, i)])[0][0][1]

            def wedge(r_in, r_out, a_c, a_h):
                cp = occ.addPoint(0, 0, 0)
                pts = [occ.addPoint(r * mm * cos(a), r * mm * sin(a), 0)
                       for a in (a_c - a_h, a_c + a_h) for r in (r_in, r_out)]
                loop = occ.addCurveLoop([occ.addLine(pts[0], pts[1]), occ.addCircleArc(pts[1], cp, pts[3]),
                                         occ.addLine(pts[3], pts[2]), occ.addCircleArc(pts[2], cp, pts[0])])
                return occ.addPlaneSurface([loop])

            # Parallel slots are boxes whose ends would be tangent to the slot-floor and
            # slot-back circles; leaving those circles in makes cusps of needle elements
            # (Newton then never converges). Both sides of them are stator iron anyway.
            radii = sorted(set(v for k, v in R.items() if not (d.slot_w and k in ("slot_a", "slot_b"))))
            surf = [annulus(0, radii[0])] + [annulus(a, b) for a, b in zip(radii[:-1], radii[1:])]
            sp_, pp_ = tau / d.slots, tau / d.poles
            tip = sorted((R["bore"], R["slot_b"] if d.outrunner else R["slot_a"]))
            for k in range(d.slots):
                if d.slot_w:
                    box = occ.addRectangle(R["slot_a"] * mm, -d.slot_w / 2 * mm, 0,
                                           (R["slot_b"] - R["slot_a"]) * mm, d.slot_w * mm)
                    lo = R["bore"] - 0.3                       # opening: a box through the bore circle
                    gate = occ.addRectangle(lo * mm, -self.open_w / 2 * mm, 0, (R["slot_a"] - lo) * mm, self.open_w * mm)
                    occ.rotate([(2, box), (2, gate)], 0, 0, 0, 0, 0, 1, k * sp_)
                    surf += [box, gate]
                else:
                    surf.append(wedge(R["slot_a"], R["slot_b"], k * sp_, 0.5 * d.slot_frac * sp_))
                    surf.append(wedge(tip[0], tip[1], k * sp_, 0.5 * self.open_frac * sp_))
            for k in range(d.poles):
                surf.append(wedge(R["mag_a"], R["mag_b"], k * pp_, 0.5 * self.magnet_arc * pp_))
            occ.synchronize()
            occ.fragment([(2, s) for s in surf], [])
            occ.synchronize()
            gmsh.model.addPhysicalGroup(2, [t for _, t in gmsh.model.getEntities(2)])
            rmid = 0.5 * sum(self.gap) * mm
            f = gmsh.model.mesh.field.add("MathEval")
            gmsh.model.mesh.field.setString(
                f, "F", f"{h_coarse*mm} - {(h_coarse-h_gap)*mm}*exp(-((sqrt(x*x+y*y)-{rmid})^2)/(2*{(2.0*mm)**2}))")
            gmsh.model.mesh.field.setAsBackgroundMesh(f)
            for opt in ("MeshSizeExtendFromBoundary", "MeshSizeFromPoints", "MeshSizeFromCurvature"):
                gmsh.option.setNumber("Mesh." + opt, 0)
            gmsh.option.setNumber("Mesh.Algorithm", 6)
            gmsh.model.mesh.generate(2)
            tags, xyz, _ = gmsh.model.mesh.getNodes()
            xy = xyz.reshape(-1, 3)[:, :2] / mm
            idmap = np.full(int(tags.max()) + 1, -1)
            idmap[tags.astype(int)] = np.arange(len(tags))
            tris = None
            for et, _, nodes in zip(*gmsh.model.mesh.getElements(2)):
                if et == 2:
                    tris = idmap[nodes.reshape(-1, 3).astype(int)]
        finally:
            gmsh.finalize()
        # weld duplicates, drop collapsed triangles and unreferenced nodes (motor/fea.py)
        uniq, inv = np.unique(np.round(xy, 5), axis=0, return_inverse=True)
        tris = inv.reshape(-1)[tris]
        tris = tris[(tris[:, 0] != tris[:, 1]) & (tris[:, 1] != tris[:, 2]) & (tris[:, 0] != tris[:, 2])]
        used = np.unique(tris)
        remap = np.full(len(uniq), -1)
        remap[used] = np.arange(len(used))
        self.xy, self.tri = uniq[used] * mm, remap[tris]          # metres

    def _materials(self):
        x, y = self.xy[self.tri, 0], self.xy[self.tri, 1]          # (M, 3)
        b = np.roll(y, -1, 1) - np.roll(y, -2, 1)
        c = np.roll(x, -2, 1) - np.roll(x, -1, 1)
        area2 = x[:, 0] * b[:, 0] + x[:, 1] * b[:, 1] + x[:, 2] * b[:, 2]
        flip = area2 < 0
        self.tri[flip] = self.tri[flip][:, [0, 2, 1]]
        if flip.any():
            return self._materials()
        self.area = area2 / 2
        self.gx, self.gy = b / area2[:, None], c / area2[:, None]  # shape-function gradients
        cx, cy = x.mean(1), y.mean(1)
        self.rc, self.th = np.hypot(cx, cy), np.arctan2(cy, cx)
        self.region = self.classify(cx * 1e3, cy * 1e3)
        self.iron = (self.region == IRON_S) | (self.region == IRON_R)
        mag = self.region == MAGNET
        pol = np.where(np.round(self.th[mag] / (tau / self.d.poles)).astype(int) % 2 == 0, 1.0, -1.0)
        self.brx, self.bry = np.zeros(len(cx)), np.zeros(len(cx))
        self.brx[mag], self.bry[mag] = self.br * pol * np.cos(self.th[mag]), self.br * pol * np.sin(self.th[mag])
        self.nu_lin = np.full(len(cx), NU0)
        self.nu_lin[mag] = NU0 / self.mur
        r_node = np.hypot(self.xy[:, 0], self.xy[:, 1])
        self.fixed = r_node > self.r["outer"] * 1e-3 - 1e-6
        self.b2_tab, self.nu_tab, self.dnu_tab, self.b_knee = _bh(self.steel)

    def _winding(self):
        """One coil per tooth. Tooth k sits between slot k and slot k+1; its coil goes
        up the near half of slot k and down the near half of slot k+1. Phase and sense
        from the tooth's electrical angle (star of slots)."""
        d = self.d
        sp_ = tau / d.slots
        ang = ((np.arange(d.slots) + 0.5) * sp_ * d.poles / 2) % tau
        sector = (np.round(ang / (pi / 3)).astype(int)) % 6
        phase = np.array([0, 2, 1, 0, 2, 1])[sector]
        sense = np.array([1, -1, 1, -1, 1, -1])[sector]
        coil = self.region == COIL
        th = self.th[coil] % tau
        k = np.round(th / sp_).astype(int) % d.slots                # slot index
        upper = ((th - k * sp_ + pi) % tau - pi) > 0                # half nearer tooth k
        tooth = np.where(upper, k, (k - 1) % d.slots)
        self.coil_idx = np.where(coil)[0]
        self.coil_phase = phase[tooth]
        self.coil_sense = sense[tooth] * np.where(upper, 1.0, -1.0)
        half = tooth * 2 + upper
        area = np.zeros(2 * d.slots)
        np.add.at(area, half, self.area[coil])
        self.coil_area = area[half]

    # ---- solve ----------------------------------------------------------------
    def solve(self, amp_turns=0.0, gamma=0.0, tol=1e-5, max_iter=40):
        """Field at `amp_turns` (peak, per coil) and electrical angle gamma. Returns
        torque (N·m) and the element flux density."""
        d = self.d
        i_ph = amp_turns * np.cos(gamma - np.array([0, tau / 3, 2 * tau / 3]))
        jz = np.zeros(len(self.area))
        jz[self.coil_idx] = self.coil_sense * i_ph[self.coil_phase] / self.coil_area
        n, tri, ar = len(self.xy), self.tri, self.area
        rhs_e = (ar * jz / 3)[:, None] + (ar * self.nu_lin)[:, None] * (self.brx[:, None] * self.gy - self.bry[:, None] * self.gx)
        rhs = np.bincount(tri.ravel(), rhs_e.ravel(), n)
        rows, cols = np.repeat(tri, 3, 1).ravel(), np.tile(tri, (1, 3)).ravel()
        free = ~self.fixed
        gg = self.gx[:, :, None] * self.gx[:, None, :] + self.gy[:, :, None] * self.gy[:, None, :]
        a = np.zeros(n)
        nu = self.nu_lin.copy()
        for it in range(max_iter):
            ae = a[tri]
            bx_, by_ = (ae * self.gx).sum(1), (ae * self.gy).sum(1)      # grad A
            b2 = bx_ * bx_ + by_ * by_
            dnu = np.zeros(len(ar))
            nu[self.iron] = np.interp(b2[self.iron], self.b2_tab, self.nu_tab)
            dnu[self.iron] = np.interp(b2[self.iron], self.b2_tab, self.dnu_tab)
            q = bx_[:, None] * self.gx + by_[:, None] * self.gy          # G^T grad A
            res = np.bincount(tri.ravel(), ((ar * nu)[:, None] * q).ravel(), n) - rhs
            ke = (ar * nu)[:, None, None] * gg + (2 * ar * dnu)[:, None, None] * q[:, :, None] * q[:, None, :]
            K = sp.csr_matrix((ke.ravel(), (rows, cols)), shape=(n, n))
            da = np.zeros(n)
            da[free] = spsolve(K[free][:, free].tocsc(), -res[free])
            step = 1.0
            r0 = np.linalg.norm(res[free])
            while step > 0.05:                                           # backtrack
                at = a + step * da
                ae = at[tri]
                tx, ty = (ae * self.gx).sum(1), (ae * self.gy).sum(1)
                nut = nu.copy()
                nut[self.iron] = np.interp((tx * tx + ty * ty)[self.iron], self.b2_tab, self.nu_tab)
                rt = np.bincount(tri.ravel(), ((ar * nut)[:, None] * (tx[:, None] * self.gx + ty[:, None] * self.gy)).ravel(), n) - rhs
                if np.linalg.norm(rt[free]) < r0 or it == 0:
                    break
                step *= 0.5
            a = a + step * da
            if np.linalg.norm(step * da) < tol * max(np.linalg.norm(a), 1e-12):
                break
        else:
            print(f"fea: Newton did not converge in {max_iter} iterations at {amp_turns:.0f} A-turns "
                  f"(check the mesh for needle elements)", file=sys.stderr)
        ae = a[tri]
        gxa, gya = (ae * self.gx).sum(1), (ae * self.gy).sum(1)
        bx, by = gya, -gxa
        cth, sth = np.cos(self.th), np.sin(self.th)
        br, bt = bx * cth + by * sth, -bx * sth + by * cth
        g = self.region == GAP
        r_i, r_o = (x * 1e-3 for x in self.gap)
        torque = d.stack * 1e-3 / (MU0 * (r_o - r_i)) * np.sum((self.rc * br * bt * ar)[g])
        return dict(torque=torque, b=np.hypot(bx, by), iters=it + 1, br_gap=br[g], th_gap=self.th[g])

    def best_gamma(self, amp_turns, n=12):
        """Electrical angle of maximum torque at this rotor position (coarse sweep,
        then a parabola through the best three)."""
        gs = np.linspace(0, tau, n, endpoint=False)
        ts = np.array([self.solve(amp_turns, g)["torque"] for g in gs])
        i = int(np.argmax(ts))
        y0, y1, y2 = ts[i - 1], ts[i], ts[(i + 1) % n]
        return gs[i] + 0.5 * (y0 - y2) / (y0 - 2 * y1 + y2) * (tau / n), ts

    def b1_gap(self):
        """Peak of the air-gap field's working harmonic at no current."""
        s = self.solve(0.0)
        pp = self.d.poles / 2
        w = self.area[self.region == GAP]
        c = np.sum(s["br_gap"] * np.cos(pp * s["th_gap"]) * w) / np.sum(w)
        s_ = np.sum(s["br_gap"] * np.sin(pp * s["th_gap"]) * w) / np.sum(w)
        return 2 * np.hypot(c, s_)

    def curve(self, amp_turns):
        """Torque at each amp-turn level, on the angle that is best at the lowest and
        re-checked at the highest."""
        g0, _ = self.best_gamma(amp_turns[0])
        g1 = max((g0 - 0.26, g0, g0 + 0.26), key=lambda g: self.solve(amp_turns[-1], g)["torque"])
        out = []
        for at in amp_turns:
            g = g0 + (g1 - g0) * (at - amp_turns[0]) / (amp_turns[-1] - amp_turns[0])
            out.append(max(self.solve(at, g)["torque"], self.solve(at, g0)["torque"]))
        return np.array(out)


def amp_turns_of(d, j):
    """Peak amp-turns per coil at copper current density j (rms)."""
    return sqrt(2) * j * d.cu_area * d.turns


if __name__ == "__main__":
    import time
    d = dz.Design()
    t0 = time.time()
    m = Motor(d)
    print(f"mesh: {len(m.xy)} nodes, {len(m.tri)} triangles, {time.time()-t0:.1f} s; iron knee at {m.b_knee:.2f} T")
    for reg, name in ((IRON_S, "stator iron"), (IRON_R, "rotor iron"), (COIL, "coil"), (GAP, "gap"), (MAGNET, "magnet")):
        print(f"  {name:12} {np.sum(m.area[m.region == reg])*1e6:8.1f} mm²")
    print(f"  slot area: mesh {np.sum(m.area[m.region == COIL])/d.slots*1e6:.2f} mm², design {d.slot_area*1e6:.2f} mm²")
    t0 = time.time()
    print(f"air-gap fundamental {m.b1_gap():.3f} T (design assumed {d.b1(96):.3f}); {time.time()-t0:.1f} s")
    at = amp_turns_of(d, 5e6)
    g, ts = m.best_gamma(at)
    print("torque vs angle at 5 A/mm²:", np.round(ts, 2))
    print(f"best angle {np.degrees(g):.0f} deg; FEA {m.solve(at, g)['torque']:.2f} N·m, "
          f"design.py unsaturated {d.kt_peak(96) * at / d.turns:.2f} N·m")
