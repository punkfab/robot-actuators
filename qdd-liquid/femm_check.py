"""
Cross-check of fea.py against FEMM 4.2 (D. Meeker, femm.info) on the shallow-slot design.

The same cross-section is drawn in FEMM from the same design.Design (same radii, slots,
magnets, winding rule, B-H points and current angle), solved with FEMM's own mesher
(Triangle) and nonlinear solver, and its weighted-stress-tensor torque on the rotor is
compared with fea.py's Arkkio torque. What is shared is the geometry and the material
data; the mesh, the solver and the torque method are independent.

FEMM is a Windows program. This runs it under wine, headless:

    WINEPREFIX=~/.wine-femm, FEMM installed at C:\\femm42 (override with FEMM_EXE)

    python qdd-liquid/femm_check.py        (or: make qdd-liquid-femm; a few minutes)
"""

import json
import os
import subprocess
import sys
from math import asin, cos, pi, sin, sqrt, tau
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import design as dz  # noqa: E402
import fea  # noqa: E402

OUT = HERE / "out" / "femm"
WINEPREFIX = os.environ.get("WINEPREFIX", str(Path.home() / ".wine-femm"))
FEMM_EXE = os.environ.get("FEMM_EXE", str(Path(WINEPREFIX) / "drive_c/femm42/bin/femm.exe"))
J_CHECK = np.array([3, 14, 26, 32.5, 43, 67, 100]) * 1e6        # A/m2 rms


def wpath(p):
    return "Z:" + str(p).replace("/", "\\\\")


def run_lua(lua, name, timeout=1500):
    """Run a Lua script in FEMM under wine with no display. A Lua error leaves FEMM
    waiting on a message box, so on timeout the hidden screen is saved for reading."""
    OUT.mkdir(parents=True, exist_ok=True)
    script = OUT / f"{name}.lua"
    script.write_text(lua)
    env = dict(os.environ, WINEPREFIX=WINEPREFIX, WINEDEBUG="-all")
    auth, disp = OUT / "xauth", 90 + os.getpid() % 9
    cmd = ["xvfb-run", "-n", str(disp), "-f", str(auth), "wine", FEMM_EXE, f"-lua-script={wpath(script)}", "-windowhide"]
    try:
        subprocess.run(cmd, env=env, check=True, capture_output=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        shot = OUT / f"{name}-stuck.png"
        subprocess.run(["import", "-display", f":{disp}", "-window", "root", str(shot)], env=dict(env, XAUTHORITY=str(auth)))
        subprocess.run(["pkill", "-x", "femm.exe"])
        raise RuntimeError(f"FEMM did not finish {script.name}; its screen is in {shot}")


def rot(u, v, a):
    return u * cos(a) - v * sin(a), u * sin(a) + v * cos(a)


class Draw:
    """Collects FEMM drawing commands for one problem type (prefix mi or hi)."""

    def __init__(self, prefix):
        self.p, self.cmd = prefix, []

    def node(self, x, y):
        self.cmd.append(f"{self.p}_addnode({x:.6f},{y:.6f})")

    def seg(self, a, b):
        for q in (a, b):
            self.node(*q)
        self.cmd.append(f"{self.p}_addsegment({a[0]:.6f},{a[1]:.6f},{b[0]:.6f},{b[1]:.6f})")

    def arc(self, r, a0, a1, maxseg=1.0):
        """Arc of radius r from angle a0 to a1 (radians, counter-clockwise, under 180 deg)."""
        a, b = (r * cos(a0), r * sin(a0)), (r * cos(a1), r * sin(a1))
        for q in (a, b):
            self.node(*q)
        self.cmd.append(f"{self.p}_addarc({a[0]:.6f},{a[1]:.6f},{b[0]:.6f},{b[1]:.6f},{np.degrees(a1 - a0):.6f},{maxseg})")

    def circle(self, r, maxseg=1.0):
        for k in range(4):
            self.arc(r, k * pi / 2, (k + 1) * pi / 2, maxseg)

    def label(self, x, y, material, mesh=0.0, group=0, magdir=None):
        self.cmd.append(f"{self.p}_addblocklabel({x:.6f},{y:.6f})")
        self.cmd.append(f"{self.p}_selectlabel({x:.6f},{y:.6f})")
        auto = 1 if mesh == 0 else 0
        if self.p == "mi":
            md = f'"{magdir}"' if isinstance(magdir, str) else (magdir or 0)
            self.cmd.append(f'mi_setblockprop("{material}",{auto},{mesh},"<None>",{md},{group},1)')
        else:
            self.cmd.append(f'hi_setblockprop("{material}",{auto},{mesh},{group})')
        self.cmd.append(f"{self.p}_clearselected()")


def bh_points(steel="silicon steel"):
    """The B-H curve fea.py uses, sampled for FEMM (which splines between points)."""
    b2, nu, _, knee = fea._bh(steel)
    b = np.concatenate([np.linspace(0.05, 1.2, 12), np.linspace(1.3, knee, 24), np.linspace(knee + 0.1, 4.0, 12)])
    h = np.interp(b * b, b2, nu) * b
    return [(0.0, 0.0)] + list(zip(b, h))


def magnetics_lua(m, cases, fem, ans_out):
    d, R = m.d, m.r
    sp_, pp_ = tau / d.slots, tau / d.poles
    g = Draw("mi")
    c = g.cmd
    c += ["newdocument(0)", f'mi_probdef(0,"millimeters","planar",1e-8,{d.stack},30)', "mi_smartmesh(1)"]
    c += ['mi_addmaterial("Air",1,1,0,0,0,0,0,1,0,0,0)', 'mi_addmaterial("Steel",1000,1000,0,0,0,0,0,1,0,0,0)']
    c += [f'mi_addbhpoint("Steel",{b:.5f},{h:.4f})' for b, h in bh_points(m.steel)]
    hc = m.br / (fea.MU0 * m.mur)
    c.append(f'mi_addmaterial("Magnet",{m.mur},{m.mur},{hc:.2f},0,0,0,0,1,0,0,0)')
    for p in range(3):
        for s in "pn":
            c.append(f'mi_addmaterial("c{p}{s}",1,1,0,0,0,0,0,1,0,0,0)')
    c.append('mi_addboundprop("zero",0,0,0,0,0,0,0,0,0)')

    # rotor: yoke ring, magnets on it
    g.circle(R["inner"])
    edges = sorted(k * pp_ + s * 0.5 * m.magnet_arc * pp_ for k in range(d.poles) for s in (-1, 1))
    for a0, a1 in zip(edges, edges[1:] + [edges[0] + tau]):
        g.arc(R["mag_a"], a0, a1)
    for k in range(d.poles):
        a0, a1 = k * pp_ - 0.5 * m.magnet_arc * pp_, k * pp_ + 0.5 * m.magnet_arc * pp_
        g.arc(R["mag_b"], a0, a1, 0.5)
        for a in (a0, a1):
            g.seg((R["mag_a"] * cos(a), R["mag_a"] * sin(a)), (R["mag_b"] * cos(a), R["mag_b"] * sin(a)))
    # stator: bore broken by the slot openings, parallel-sided slots split into two coil sides
    ow, sw, rb, sa, sb = m.open_w, d.slot_w, R["bore"], R["slot_a"], R["slot_b"]
    da, ub = asin(ow / 2 / rb), sqrt(rb * rb - ow * ow / 4)
    for k in range(d.slots):
        a = k * sp_
        g.arc(rb, a + da, a + sp_ - da, 0.5)
        P = lambda u, v: rot(u, v, a)
        outline = [(ub, ow / 2), (sa, ow / 2), (sa, sw / 2), (sb, sw / 2), (sb, 0), (sb, -sw / 2), (sa, -sw / 2),
                   (sa, -ow / 2), (ub, -ow / 2)]
        for q0, q1 in zip(outline, outline[1:]):
            g.seg(P(*q0), P(*q1))
        g.seg(P(sa, -ow / 2), P(sa, 0))
        g.seg(P(sa, 0), P(sa, ow / 2))
        g.seg(P(sa, 0), P(sb, 0))
    g.circle(R["outer"], 2.0)
    for k in range(4):
        a = (k + 0.5) * pi / 2
        c.append(f"mi_selectarcsegment({R['outer']*cos(a):.5f},{R['outer']*sin(a):.5f})")
    c += ['mi_setarcsegmentprop(2,"zero",0,0)', "mi_clearselected()"]

    g.label(0, 0, "Air", 3.0)
    g.label((R["inner"] + R["mag_a"]) / 2, 0, "Steel", fem["iron"], group=1)
    for k in range(d.poles):
        r = (R["mag_a"] + R["mag_b"]) / 2
        g.label(r * cos(k * pp_), r * sin(k * pp_), "Magnet", fem["magnet"], group=1,
                magdir="theta" if k % 2 == 0 else "theta+180")
    rg = sum(m.gap) / 2
    g.label(rg * cos(sp_ / 2), rg * sin(sp_ / 2), "Air", fem["gap"])
    ry = (sb + R["outer"]) / 2
    g.label(ry * cos(sp_ / 2), ry * sin(sp_ / 2), "Steel", fem["iron"])
    # winding: the rule in fea.Motor._winding
    ang = ((np.arange(d.slots) + 0.5) * sp_ * d.poles / 2) % tau
    sector = (np.round(ang / (pi / 3)).astype(int)) % 6
    phase, sense = np.array([0, 2, 1, 0, 2, 1])[sector], np.array([1, -1, 1, -1, 1, -1])[sector]
    for k in range(d.slots):
        for upper in (True, False):
            tooth = k if upper else (k - 1) % d.slots
            s = sense[tooth] * (1 if upper else -1)
            x, y = rot((sa + sb) / 2, (sw / 4) * (1 if upper else -1), k * sp_)
            g.label(x, y, f"c{phase[tooth]}{'p' if s > 0 else 'n'}", fem["coil"])
    c.append(f'mi_saveas("{wpath(OUT / "motor.fem")}")')
    c.append(f'f = openfile("{wpath(ans_out)}","w")')
    half = (sw / 2) * (sb - sa) * 1e-6                              # m2, one coil side
    for name, at, gamma in cases:
        i_ph = at * np.cos(gamma - np.array([0, tau / 3, 2 * tau / 3]))
        for p in range(3):
            j = i_ph[p] / half / 1e6                                 # MA/m2
            c += [f'mi_modifymaterial("c{p}p",4,{j:.6f})', f'mi_modifymaterial("c{p}n",4,{-j:.6f})']
        c += ["mi_analyze(1)", "mi_loadsolution()", "mo_groupselectblock(1)", "t = mo_blockintegral(22)",
              "mo_clearblock()", f'write(f,"T {name} ",t,"\\n")']
        if at == 0:                                                  # air-gap field along the mid-gap circle
            for a in np.arange(1440) * tau / 1440:               # (Lua 4 trig is in degrees: do it here)
                x, y = rg * cos(a), rg * sin(a)
                c.append(f'A,bx,by = mo_getpointvalues({x:.6f},{y:.6f}); write(f,"B {a:.6f} ",bx*{cos(a):.8f}+by*{sin(a):.8f},"\\n")')
        c += ['write(f,"N ",mo_numelements(),"\\n")', "mo_close()"]
    c += ["closefile(f)", "quit()"]
    return "\n".join(c) + "\n"


def run_magnetics(d=dz.SHALLOW, fem=None, tag="motor"):
    fem = fem or dict(gap=0.25, iron=1.2, magnet=0.8, coil=1.2)
    m = fea.Motor(d)
    ats = [fea.amp_turns_of(d, j) for j in J_CHECK]
    g0, _ = m.best_gamma(ats[0])
    cases = [("open", 0.0, 0.0)] + [(f"j{k}", at, g0) for k, at in enumerate(ats)]
    cases += [(f"g{k}", ats[3], g0 + dg) for k, dg in enumerate((-0.5, -0.25, 0.25, 0.5))]
    ans = OUT / f"{tag}.txt"
    run_lua(magnetics_lua(m, cases, fem, ans), tag)
    t, b, n = {}, [], 0
    for line in ans.read_text().split("\n"):
        w = line.split()
        if w[:1] == ["T"]:
            t[w[1]] = float(w[2])
        elif w[:1] == ["B"]:
            b.append((float(w[1]), float(w[2])))
        elif w[:1] == ["N"]:
            n = int(float(w[1]))
    b = np.array(b)
    pp = d.poles / 2
    b1 = 2 * np.hypot(np.mean(b[:, 1] * np.cos(pp * b[:, 0])), np.mean(b[:, 1] * np.sin(pp * b[:, 0])))
    rows = [dict(case=name, amp_turns=at, gamma=g, ours=float(m.solve(at, g)["torque"]), femm=t[name])
            for name, at, g in cases]
    return dict(rows=rows, b1_ours=float(m.b1_gap()), b1_femm=float(b1), tri_ours=int(len(m.tri)), tri_femm=n, gamma0=float(g0))


def main():
    r = run_magnetics()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "magnetics.json").write_text(json.dumps(r, indent=1))
    print(f"mesh: fea.py {r['tri_ours']} triangles, FEMM {r['tri_femm']}")
    print(f"no-load air-gap fundamental: fea.py {r['b1_ours']:.4f} T, FEMM {r['b1_femm']:.4f} T "
          f"({(r['b1_ours']/r['b1_femm']-1)*100:+.2f}%)")
    print(f"\n  {'case':28} {'fea.py N·m':>11} {'FEMM N·m':>10} {'difference':>11}")
    names = ["no current (cogging at this position)"] + [f"{j/1e6:g} A/mm²" for j in J_CHECK] \
        + [f"32.5 A/mm², angle {dg:+.2f} rad" for dg in (-0.5, -0.25, 0.25, 0.5)]
    for name, row in zip(names, r["rows"]):
        diff = f"{(row['ours']/row['femm']-1)*100:+.2f}%" if abs(row["femm"]) > 0.5 else f"{row['ours']-row['femm']:+.3f} N·m"
        print(f"  {name:28} {row['ours']:11.3f} {row['femm']:10.3f} {diff:>11}")


if __name__ == "__main__":
    main()
