# RobStride: designing with QDD actuators

RobStride's RS series (RS00–RS06, RS10P) are quasi-direct-drive (QDD) joints: an
outrunner BLDC, one planetary stage (7.75–10:1, 25:1 for the RS10P) and a FOC driver on
CAN. They descend from Ben Katz's mini-cheetah actuator. Their "MIT mode" is Katz's frame:
16-bit p* and 12-bit v*, kp, kd, τff, running the law

    τ = kp·(p* − q) + kd·(v* − q̇) + τff          on the drive, per current-loop cycle

This folder turns them into something a design can use: data, a MuJoCo model, a sizing
tool, CAD, and the bench loop to replace assumptions with measurements.

| file | what |
|---|---|
| `specs.py` | every model's numbers, each with its source (booklet, manuals, sysid, Katz) |
| `thermal.py` | two-node winding model fitted to the booklet's overload tables |
| `servo.py` | `motorize()` + `RobStrideServos`: the MIT law in MuJoCo, with host ZOH, latency, packet quantization, the T-N envelope, gearbox friction, rotor armature, thermal trip |
| `sizing.py` | lightest RobStride per joint from inverse dynamics: peak, thermal, T-N; actuator masses fed back |
| `demo.py` | size → tune → verify on a 3-DOF arm with 2 kg payload (`make robstride-demo`) |
| `bench.py` / `record.py` | pendulum-bench sysid: record over CAN, regress, polish, validate |
| `can_mit.py` | MIT frames with each model's own field ranges |
| `cad.py` | the vendor STEPs, normalized: joint axis +Z, output flange face at z = 0 |

## Sources

- **RobStride's own GitHub** ([Product_Information](https://github.com/RobStride/Product_Information)):
  - the spec booklet: T-N tables, overload tables, R/L, output-side inertia, drawings
  - per-model manuals: MIT field ranges
  - STEP files
- **[Robstride-SysID](https://github.com/T-K-233/Robstride-SysID)**: measured drag, damping and armature.
- **B. Katz, *A Low Cost Modular Actuator for Dynamic Robots*** (MIT MS thesis, 2018): the
  ancestor design, with measured gearbox friction, backlash, loop rates and the frame.
- **Wensing et al., *Proprioceptive Actuator Design in the MIT Cheetah*** (T-RO 2017): the
  design rationale (see [Notes from the papers](#notes-from-the-papers)).

Anything assumed is marked `[ASSUMED]` in `specs.py`.

## Using it in a design

```python
from servo import motorize, RobStrideServos
from sizing import size_joints

picks, checks, every = size_joints(make_spec, joints, qpos_traj, dt)   # {"hip": "RS03", ...}
motorize(spec, picks)                                   # joints -> torque motors
m = spec.compile(); d = mujoco.MjData(m)
rs = RobStrideServos(m, d, picks, kp=80, kd=4, host_hz=500)
rs.set("hip", p=0.4, v=0.0, tau_ff=tau_id)             # every host tick
rs.step()                                              # law + friction + mj_step
rs.state("hip")                                        # q, dq, tau, theta (1 = trip)
```

For assemblies, `cad.actuator("RS03", loc)` gives the vendor model placed with the output
flange on your link; `cad.envelope()` is the no-network cylinder.

## What the demo shows (`make robstride-demo`)

The design is a 3-DOF arm (0.35 + 0.30 + 0.10 m) carrying 2 kg through a pick-and-place
cycle, 3.2 s including 0.4 s dwells.

**1. Heat decides the size, not peak torque.** The shoulder peaks at 45 N·m, and an RS03's
60 N·m covers that. But the cycle's RMS is 21 N·m, which is right at the RS03's rating, and
the dwells count against the lower *holding* rating (13 N·m: one phase carries the current
when the rotor is still). So the pick is an RS04 at 39% thermal load:

- shoulder RS04
- elbow RS10P
- wrist RS00
- 2.19 kg of actuators, fed back into the inverse dynamics

**2. A datasheet guess fails.** Choosing "peak ≥ 1.5× static gravity" and ignoring the
actuators' own mass picks RS10P / RS00 / RS05. In the model, the elbow's RS00 hits its
over-temperature trip after **8 s**, the driver cuts torque, and the arm drops.

**3. Tracking comes from feed-forward, and what's left is friction and host rate.** Tool-tip
error over the last cycle:

| | mean | max |
|---|---|---|
| ideal MuJoCo servo, PD only | 16.5 mm | 50 mm |
| RobStride model, PD only | 16.6 mm | 54 mm |
| ideal servo, PD + inverse-dynamics FF | 0.4 mm | 1.1 mm |
| **RobStride model, PD + FF** | **1.8 mm** | **4.5 mm** |
| same, randomized (friction ×0.5–2, latency 0.3–2 ms, host 250–1000 Hz, 40–50 V, armature ×1–3) | | 4.6–10.9 mm |

Unlike the STS3215 (`../sts3215`), the MIT law hides nothing: kp *is* the stiffness. What
the ideal servo leaves out is the gearbox friction, which the feed-forward can't know, and
the host's zero-order hold. The ablation behind this:

- With zero latency and a 1 kHz host, friction alone takes the mean error from 0.18 to
  1.0 mm.
- At a tenth of the gains, friction takes it from 1.2 to 18 mm.
- Doubling friction takes the max from 4.5 to 9.1 mm.
- A 100 Hz host gives 5.2 mm.

**4. The host rate shows up as torque ripple.** The drive holds p* between frames (Katz's
firmware doesn't interpolate), so each frame steps the torque by kp·v/f_host. On the
shoulder, at kp 2100 and 500 Hz, that's about 10 N·m of sawtooth during fast moves. This is
worth checking on hardware; it sets a floor on host rate for stiff gains.

(A trap I hit building this: computing the feed-forward with `mj_inverse` on the *simulated*
model includes the friction the updater just set, so the "controller" knew the true
friction every step. The host model has to be a separate copy without it.)

## Before hardware: two things to know

**lerobot's RobStride driver and the manuals disagree on field ranges.** The MIT fields are
scaled per model by the firmware:

| model | kp | kd | v |
|---|---|---|---|
| RS00/01/02/05/10P | 0–500 | 0–5 | model-specific |
| RS03/04/06 | 0–5000 | 0–100 | RS03 ±20, RS04 ±15 |

lerobot's `motors/robstride` (2026-06) encodes kp over 0–500 and kd over 0–5 for every
model, and v over ±33 for RS02/03/04. If the manuals are right, lerobot's kp 40 / kd 2 would
arrive at an RS03 as **kp 400 / kd 40** (`python can_mit.py` prints this per model). Not
verified on hardware here. `can_mit.encode()` uses the manual ranges.

**What RobStride doesn't publish:**
- gearbox efficiency
- backlash
- loop rates

The prior uses Katz's measured form for the ancestor design:
- drag + 0.04·|τ_motor|, so η ≈ 0.96 both ways
- drag per model from Robstride-SysID
- backlash about 0.005 rad (not modelled yet)

The demo's band shows that the friction assumption is the biggest open variable, so measure it:

## The bench (`make robstride-bench`)

**The rig.**
- A 0.25 m arm on the output, axis horizontal.
- Masses giving 25%, 40% and 60% of rated torque. Two loads are needed: at one, drag and
  load-proportional friction are nearly interchangeable.

**The trajectories**, shared between `record.py` (hardware) and the sim:

| trajectory | what it does |
|---|---|
| coast | release from 90° and let it swing down |
| sweep | torque chirp |
| raise | slow position triangle |
| steps | position steps (held out for validation) |

**The fit.**
1. An equation-error regression. Armature, drag, Km, Ke and Kv all enter the pendulum's
   dynamics linearly.
2. Nelder–Mead on the replayed trajectories, to polish.

Fitting trajectories from the prior alone stalled at 0.8°, because long swings make a bumpy
objective.

**Self-test on RS03.** The "real" unit has drag ×1.6, η 92%/92.5% and armature ×2.2:

| | drag | η drive / back | Kv | armature |
|---|---|---|---|---|
| truth | 0.320 | 92.2% / 92.5% | 0.0200 | 0.0440 |
| regression | 0.336 | 92.1% / 92.4% | 0.0174 | 0.0469 |
| **fitted** | **0.320** | **92.1% / 92.4%** | **0.0200** | **0.0440** |

Training MAE goes 8.7° → 0.01°, and the held-out steps 0.21° → 0.01°.
Pass the result to `RobStrideServos(friction=..., armature=...)`.

**On hardware** (`record.py`, not yet run: `--check` first):

```bash
pip install python-can
sudo ip link set can0 up type can bitrate 1000000       # motor in MIT protocol mode
python robstride/record.py --kind RS03 --id 1 --check
python robstride/record.py --kind RS03 --id 1 --zero    # arm hanging
python robstride/record.py --kind RS03 --id 1 --load 0.25 --traj coast   # ... the bench.PLAN set
python robstride/bench.py RS03 ~/robstride_logs/*.npz
```

## CAD (`make robstride-cad`)

- `cad.py` fetches the STEPs into `robstride/cad/`. It doesn't vendor them: RobStride states no
  licence.
- It finds each model's joint axis, and finds the output end from the tapped holes on the
  booklet's output bolt circle. The vendor files disagree: four have the output on +Z, four on −Z.
- It places the flange face at z = 0 and caches the result as BREP. STEP import takes about a
  minute per model; the cached BREP loads in under a second.
- All eight match the booklet's OD×H.

## Notes from the papers

- **Katz 2018** (the ancestor):
  - The PD law runs on the drive at the 40 kHz current-loop rate. A CAN timeout zeroes the
    current, so the motor goes limp.
  - The torque estimate is a constant Kt·iq, while the real Kt drops 12% at peak current. So
    reported torque reads high near the peak.
  - Braking (negative work) gets full torque at any speed, which `envelope()` models.
  - Winding time constant is about 39 s, with a slower tail, which is why the thermal model
    has two nodes.
  - His impact rule, τ_gear ≈ ω·J_reflected·ω_n, bounds the gear load in a collision. It
    isn't a sizing check here yet.
- **Wensing et al. 2017** (MIT Cheetah):
  - Reflected rotor inertia, N²·J_r, is MuJoCo's `armature`, which is how it's modelled here.
  - The impact mitigation factor (IMF, from the operational-space inertia with and without
    it) and the effective-mass bandwidth √(k/m_e) are the design metrics for legs. Not
    implemented yet.

## Caveats

- **Nothing has touched a RobStride yet.** The self-tests show the tools recover a known truth.
  They don't show the priors are right.
- **The thermal model is fitted to vendor tables measured on aluminium heatsinks**
  (200×200 mm for RS03, 220×200 for RS04). A motor buried in a printed link runs hotter.
- **No backlash, no Kt saturation, no cogging** (Katz measured 0.14 N·m RMS at the output).
  No field weakening.
- **Some published numbers conflict** (rated torques, RS10P Kt, RS00 revision). `specs.py` notes
  each conflict and uses the conservative value.
