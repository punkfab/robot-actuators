"""
Scenarios for sim.py --demo (see sim.py's docstring). Writes out/sim_demo.png, out/sim.json and
out/sim.gif (a cutaway, slow motion).
"""

import json
import os
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
os.environ.setdefault("MUJOCO_GL", "egl")
import mujoco  # noqa: E402

import design as D  # noqa: E402
import sim as S  # noqa: E402

OUT = HERE / "out"


def single_block(w):
    return replace(w, blocks=((0.0, 0.0, False),))


def light_arm(p):
    return replace(p, arm_mass=0.01, arm_len=0.05)


# ---------------------------------------------------------------- scenarios

def ratio(dphi=100.0, T_end=0.065):
    """spin at constant voltage from the conjugate state: output/input, input torque ripple"""
    p = light_arm(S.Params())
    sim = S.Wobble(p=p)
    sim.reset(0.0, dphi)
    sim.mode, sim.volts = "volts", p.Kt * dphi
    t, phi, out, tau = [], [], [], []
    a0 = 0.0
    for k in range(int(T_end / p.dt)):
        info = sim.step()
        if k % 50 == 0:
            t.append(sim.d.time)
            phi.append(info["phi"])
            a = sim.output_angle()
            out.append(a)
            tau.append(info["tau_in"])
    phi, out = np.array(phi), np.unwrap(np.array(out))
    slope = np.polyfit(phi, out, 1)[0]
    return dict(t=np.array(t), phi=phi, out=out, tau=np.array(tau), ratio=1 / slope, expected=sim.w.ratio)


def hold_case(w=None, p=None, T_end=0.25):
    sim, info, tr = S.hold(w, p, T_end=T_end, every=100)
    F, Mo = S.arm_load(sim)
    return sim, info, tr, F, Mo


def eff(mu, dphi=150.0, T_end=0.1):
    """lift the arm at speed: efficiency from the energy balance"""
    p = replace(S.Params(), mu_roll=mu)
    sim = S.Wobble(p=p)
    sim.reset(0.0, dphi)
    T_load = 4.2 * 9.81 * 0.40
    sim.mode = "volts"
    sim.volts = p.Kt * dphi + p.R * (T_load / sim.w.ratio / 0.9 / p.Kt)
    m, d = sim.m, sim.d
    arm = m.body("arm").id
    out = sim.bid[-1]
    n_settle = int(0.02 / p.dt)
    W_m, E0 = 0.0, None
    for k in range(int(T_end / p.dt)):
        info = sim.step()
        if k == n_settle:
            E0 = energies(sim)
        if k > n_settle:
            W_m += info["tau_m"] * info["dphi"] * p.dt
    E1 = energies(sim)
    dPE_out = E1["pe_out"] - E0["pe_out"]
    dKE_out = E1["ke_out"] - E0["ke_out"]
    dE_rest = (E1["total"] - E1["pe_out"] - E1["ke_out"]) - (E0["total"] - E0["pe_out"] - E0["ke_out"])
    useful = dPE_out + dKE_out
    eta = useful / (W_m - dE_rest)
    return dict(mu=mu, eta=eta, W_m=W_m, useful=useful)


def energies(sim):
    m, d = sim.m, sim.d
    mujoco.mj_energyPos(m, d)
    mujoco.mj_energyVel(m, d)
    tot = d.energy.sum()
    ids = [b for b in range(m.nbody) if b == sim.bid[-1] or m.body_parentid[b] == sim.bid[-1]]
    pe = sum(m.body_mass[b] * 9.81 * d.xipos[b][2] for b in ids)
    ke = 0.0
    for b in ids:
        v6 = np.zeros(6)
        mujoco.mj_objectVelocity(m, d, mujoco.mjtObj.mjOBJ_BODY, b, v6, 0)
        # velocity at the body frame origin; shift to the COM
        om, vo = v6[:3], v6[3:]
        vc = vo + np.cross(om, d.xipos[b] - d.xpos[b])
        R = d.ximat[b].reshape(3, 3)
        I = R @ np.diag(m.body_inertia[b]) @ R.T
        ke += 0.5 * m.body_mass[b] * vc @ vc + 0.5 * om @ I @ om
    return dict(total=tot, pe_out=pe, ke_out=ke)


def shake(w, p, dphi=314.0, T_end=0.05):
    """3000 rpm, no load: the frame's reaction moment (CAD x, y components)"""
    p = light_arm(p)
    sim = S.Wobble(w, p)
    sim.reset(0.0, dphi)
    sim.mode, sim.volts = "volts", p.Kt * dphi
    t, Mxy, tau = [], [], []
    for k in range(int(T_end / p.dt)):
        info = sim.step()
        if k % 20 == 0:
            t.append(sim.d.time)
            Mxy.append(np.hypot(*sim.ground[3:5]))
            tau.append(info["tau_in"])
        if not np.isfinite(sim.d.qpos).all():
            break
    t, Mxy = np.array(t), np.array(Mxy)
    s = t > 0.01
    return dict(t=t, M=Mxy, rms=float(np.sqrt(np.mean(Mxy[s] ** 2))), tau=np.array(tau))


# ---------------------------------------------------------------- gif

def gif(path, revs=3, frames_per_rev=30, dphi=125.0):
    from PIL import Image, ImageDraw, ImageFont
    p = light_arm(S.Params())
    S.CUTAWAY = True
    sim = S.Wobble(p=p)
    S.CUTAWAY = False
    sim.reset(0.0, dphi)
    sim.mode, sim.volts = "volts", p.Kt * dphi
    r = mujoco.Renderer(sim.m, 540, 960)
    cam = mujoco.MjvCamera()
    cam.lookat[:] = (0.0, 0.0, 0.0)
    cam.distance, cam.azimuth, cam.elevation = 0.24, 28.0, -22.0
    per = int(2 * np.pi / frames_per_rev / dphi / p.dt)
    font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf", 20)
    ims = []
    for f in range(revs * frames_per_rev):
        sim.step(per)
        r.update_scene(sim.d, cam)
        im = Image.fromarray(r.render())
        dr = ImageDraw.Draw(im)
        dr.text((14, 10), f"twin wobble reducer {abs(sim.w.ratio):.0f}:1 - MuJoCo 4.0, roller contact law, slow motion",
                fill=(20, 20, 25), font=font)
        dr.text((14, 36), f"input {sim.d.qpos[sim.qin] / 2 / np.pi:5.2f} rev   output {np.degrees(sim.output_angle()):6.2f} deg",
                fill=(60, 60, 70), font=font)
        ims.append(im)
    ims[0].save(path, save_all=True, append_images=ims[1:], duration=50, loop=0)
    return path


# ---------------------------------------------------------------- main

def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    w = D.Wobble()
    res = {}

    print("ratio ...")
    rt = ratio()
    res["ratio"] = dict(measured=rt["ratio"], expected=rt["expected"],
                        tau_in_ripple=float(rt["tau"].max() - rt["tau"].min()))
    print(f"  output/input measured 1/{rt['ratio']:.3f}, expected 1/{rt['expected']:.3f}")

    print("hold (twin, 8 kg at 0.4 m, 0.1 m off the plane) ...")
    p_hold = replace(S.Params(), arm_mass=8.0, arm_offset=0.10)
    sim, info, tr, F, Mo = hold_case(p=p_hold)
    s = S.solve_general(w, info["phi"], F, Mo * 1e3, sim.p.k_in * 1e3)
    R = S.G.T @ sim.d.xmat[sim.bid[-1]].reshape(3, 3)
    res["hold"] = dict(load_N=F.tolist(), moment_Nm=Mo.tolist(), roller_Q_max=float(info["Q"].max()),
                       solve_Q_max=float(s["Q"].max()), tilt_arcmin=float(tr["tilt"][-1]),
                       solve_tilt_arcmin=float(np.degrees(np.linalg.norm(s["xo"][3:5])) * 60),
                       p_max_MPa=float(D.hertz(w, info["Q"], info["k"]).max()), motor_A=float(sim.i_motor))
    print(f"  load: torque {Mo[2]:.1f} N m, tilt moment {np.hypot(Mo[0], Mo[1]):.1f} N m, radial {np.hypot(F[0], F[1]):.0f} N")
    print(f"  roller max {info['Q'].max():.0f} N (design.solve {s['Q'].max():.0f}); output tilt {tr['tilt'][-1]:.2f} arcmin "
          f"(design.solve {res['hold']['solve_tilt_arcmin']:.2f}); peak pressure {res['hold']['p_max_MPa']:.0f} MPa")

    print("A/B: one block, no output bearing / one block + aux bearing / twin ...")
    p_ab = replace(S.Params(), arm_mass=4.0)
    ab = {}
    for name, ww, pp in (("twin", w, p_ab), ("one block, no bearing", single_block(w), p_ab),
                         ("one block + crossed roller", single_block(w),
                          replace(p_ab, k_aux=1e8, k_aux_tilt=D.HPF25["K_tilt"]))):
        sim_ab, info_ab, tr_ab = S.hold(ww, pp, T_end=0.15, every=100)
        ab[name] = tr_ab
        print(f"  {name:26s} output {tr_ab['ang'][-1]:+8.2f} deg, tilt {tr_ab['tilt'][-1]:8.2f} arcmin, "
              f"max roller {tr_ab['Qmax'][-1]:7.0f} N, t {tr_ab['t'][-1]:.3f} s")
    res["ab"] = {k: dict(final_deg=float(v["ang"][-1]), tilt_arcmin=float(v["tilt"][-1]), t_end=float(v["t"][-1]))
                 for k, v in ab.items()}

    print("efficiency ...")
    ef = [eff(0.0), eff(S.Params().mu_roll)]
    res["efficiency"] = ef
    for e in ef:
        print(f"  mu {e['mu']:.3f}: eta {e['eta'] * 100:.1f} %")

    print("shake at 3000 rpm ...")
    sh_twin = shake(w, S.Params())
    sh_one = shake(single_block(w), replace(S.Params(), k_aux=1e8, k_aux_tilt=D.HPF25["K_tilt"]))
    res["shake"] = dict(twin_rms_Nm=sh_twin["rms"], one_block_rms_Nm=sh_one["rms"])
    print(f"  frame reaction moment rms: twin {sh_twin['rms']:.3f} N m, one block {sh_one['rms']:.3f} N m")

    fig, ax = plt.subplots(2, 3, figsize=(16, 8.5))
    a = ax[0, 0]
    a.plot(rt["phi"] / 2 / np.pi, np.degrees(rt["out"]), label="output")
    a.plot(rt["phi"] / 2 / np.pi, np.degrees(rt["phi"] / rt["expected"]), "--", label=f"input / {rt['expected']:.0f}")
    a.set_xlabel("input (rev)")
    a.set_ylabel("output (deg)")
    a.set_title(f"ratio: measured {rt['ratio']:.2f}, expected {rt['expected']:.2f}")
    a.legend()
    a = ax[0, 1]
    a.plot(tr["t"] * 1e3, tr["tilt"], label="output tilt (arcmin)")
    a.plot(tr["t"] * 1e3, np.array(tr["Qmax"]) / 100, label="max roller load (100 N)")
    a.axhline(res["hold"]["solve_tilt_arcmin"], ls="--", c="k", lw=0.8, label="design.solve tilt")
    a.set_xlabel("ms")
    a.set_title(f"hold 8 kg at 0.4 m: T {Mo[2]:.0f} N m, M {np.hypot(Mo[0], Mo[1]):.0f} N m, no output bearing")
    a.legend(fontsize=8)
    a = ax[0, 2]
    for k, v in ab.items():
        a.plot(v["t"] * 1e3, v["ang"], label=k)
    a.set_xlabel("ms")
    a.set_ylabel("output angle (deg)")
    a.set_title("A/B: the one-block output wobbles free")
    a.legend(fontsize=8)
    a = ax[1, 0]
    a.bar([f"mu {e['mu']}" for e in ef], [e["eta"] * 100 for e in ef], color=["#999", "#d9480f"])
    a.set_ylim(0, 105)
    a.set_ylabel("efficiency (%)")
    a.set_title("lifting 16.5 N m at 150 rad/s input (preload 150 N)")
    a = ax[1, 1]
    a.plot(sh_one["t"] * 1e3, sh_one["M"], label=f"one block + crossed roller ({sh_one['rms']:.2f} rms)")
    a.plot(sh_twin["t"] * 1e3, sh_twin["M"], label=f"twin ({sh_twin['rms']:.2f} rms)")
    a.set_xlabel("ms")
    a.set_ylabel("frame moment |Mxy| (N m)")
    a.set_title("shaking moment at 3000 rpm, no load")
    a.legend(fontsize=8)
    a = ax[1, 2]
    a.plot(rt["phi"] / 2 / np.pi, rt["tau"] * 1e3)
    a.set_xlabel("input (rev)")
    a.set_ylabel("input torque from the blocks (mN m)")
    a.set_title("input torque, no load (inertia + friction)")
    fig.tight_layout()
    fig.savefig(OUT / "sim_demo.png", dpi=110)
    (OUT / "sim.json").write_text(json.dumps(res, indent=1, default=float))
    print("gif ...")
    gif(OUT / "sim.gif")
    print("wrote out/sim_demo.png, out/sim.json, out/sim.gif")


if __name__ == "__main__":
    main()
