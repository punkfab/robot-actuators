# Real2sim across punkfab

Goal: a punkfab sim's claim ("reaches the datum", "holds the load", "the finger closes")
should hold on the real hardware. Most of the gap at our scale is the **actuator**, so
that is the first thing to model honestly.

Two external references set the recipe:
- **Rhoban's [BAM](https://github.com/Rhoban/bam)** (ICRA 2025): extended friction models plus the
  servo's firmware law, identified from pendulum logs.
- **Pollen's [Microduck](https://github.com/pollen-robotics/microduck_rl)**, which shipped sim2real with
  BAM M6 servos, domain randomization over friction, voltage, sag and delay, series backlash
  that the encoder reads through, and a one-servo sim-vs-real testbench.

## The recipe, and where punkfab stands

| Step | What | Status |
|---|---|---|
| 1 | Model the actuator's control law, not an ideal servo | ✅ STS3215 via BAM (`sts3215/`); custom gearboxes via `friction/` |
| 2 | Load-dependent / directional / Stribeck friction, per step | ✅ `friction/` (matches BAM's code to 1e-15); MuJoCo creep fixed |
| 3 | Measure your own units: record → replay → fit | 🟡 pipeline built and verified on synthetic logs; **needs a real SO-101 run** |
| 4 | Randomize what you can't pin down; report a band | 🟡 friction/voltage bands in `sts3215/compare.py` and `crosscheck.py`; no delay or sag DR yet |
| 5 | Backlash in series, encoder through the play | ⬜ upstream MJCF has an unused ±0.5° `backlash` class |
| 6 | Re-run every gate that steps a real actuator with the honest model | 🟡 software-mfg workcell gate done (still passes, 4–6 mm vs 0.4 mm) |

## By project

- **SO-101 arm** (software-mfg, autobot, so101-lab, physical-ml, robot-effectors):
  - Use `bamify` + `StsServos`.
  - Next steps:
    - record real logs (`sts3215/record.py`) and fit
    - add backlash
    - re-check the rest of software-mfg's arm gates: the toolchange carry and the cell hand-off
  - autobot scores reach kinematically, so sag isn't in its ranking yet. A dynamic settle per
    candidate would add it.
- **Our own drives** (cycloidal, planetary, capstan):
  - `Friction.from_efficiency` gives a prior before printing.
  - A printable pendulum bench is the way to identify the real part. It uses the same four
    trajectory types as BAM, and BAM's fitter works on its logs directly.
- **MG996R** (tendril), **TT motors**, **steppers** (wirebender, NEMA-17 drives): no fitted
  model exists.
  - The hobby servos can go on the same bench as a new BAM actuator class.
  - Steppers need a different model (pull-out torque vs speed, missed steps), not friction.
- **Feetech wheel mode** (robot-locomotion omni base): BAM's STS3215 fit is for position
  mode. Velocity mode runs a different firmware law and needs its own log.
