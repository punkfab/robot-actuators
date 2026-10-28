# Robot actuators

Actuator designs for robot arms and legs: gear reducers, linear actuators, motors and
soft or compliant mechanisms. Each one is defined in code, simulated, and drawn as an
assembly before anything is built.

The approach: change a parameter, re-run, and see the geometry and the predicted behaviour
together. That covers torque, speed, efficiency, life and fit, fast enough to compare
designs without a dyno. Most of what's here has not been built. Each project's README says
what was checked and how.

![Exploded view of the printable RV reducer, from rv/](rv/out/explode_printed.png)

## Featured

| Project | What | Write-up |
|---|---|---|
| [`rv/`](rv/) | An RV reducer from a crankshaft-fatigue paper, built two ways. One is steel with only catalogue bearings (800 N·m). The other is 3D-printable, 81:1, with only ball bearings and dowel pins. The paper's fatigue model is rebuilt; there's an exploded drawing and an assembly animation. | [blog](https://punkfab.com/blog/rv-reducer/) |
| [`rollerscrew/`](rollerscrew/) | An Optimus-style linear actuator: a hollow motor turning a planetary roller screw, 8 kN class, sized against a day of knee duty. Plus a MuJoCo leg with the actuator in its knee, for landings and squats. | |
| [`qdd-liquid/`](qdd-liquid/) | A joint motor cooled through hollow copper conductors. Sizing, a nonlinear field solution checked against FEMM, and CAD. | [blog](https://punkfab.com/blog/liquid-cooled-qdd/), [preprint](https://doi.org/10.5281/zenodo.23155222) |
| [`friction/`](friction/), [`sts3215/`](sts3215/) | Honest servo models for simulation: extended gearbox friction (after Rhoban's BAM), and the SO-101's STS3215 servo fitted from logs. The plan across projects is in [`REAL2SIM.md`](REAL2SIM.md). | [blog](https://punkfab.com/blog/honest-actuators/) |
| [`muscle/`](muscle/) | A matrix muscle of SMA cells, with a HASEL pouch cell and a slow-manifold controller tried on both. | [blog](https://punkfab.com/blog/sma-vs-hasel/) |
| [`chain-motor/`](chain-motor/) | A two-cell Variable Chain Motor (geared BLDC cells that bend across a joint), built from the IROS 2025 paper. | [blog](https://punkfab.com/blog/paper-to-assembly/) |

## All projects

Run `make help` for every target. The main one for each project:

**Gear reducers**

| Folder | What | `make` |
|---|---|---|
| [`cycloidal/`](cycloidal/) | Hybrid planetary + cycloidal, 40:1, for a 2204 drone motor; the first build, with the efficiency model and MuJoCo bench | `cycloidal`, `actuator-bench` |
| [`cycloidal-center/`](cycloidal-center/) | 10:1 cycloidal with a straddled carrier for moment load | `cycloidal-center` |
| [`cycloidal-2stage/`](cycloidal-2stage/) | Two-stage compound cycloidal, 121:1, body output (Blagojevic) | `cycloidal-2stage` |
| [`planetary/`](planetary/) | Standalone 4:1 planetary for a 2204 motor | `planetary` |
| [`planetary-inverted/`](planetary-inverted/) | 5:1 carrier-out planetary, NEMA 17 in and out | `planetary-inverted` |
| [`rv/`](rv/) | RV reducer: fatigue model, steel and printable versions | `rv`, `rv-printed` |

**Linear actuators**

| Folder | What | `make` |
|---|---|---|
| [`rollerscrew/`](rollerscrew/) | Roller-screw actuator, 8 kN knee class, with a MuJoCo leg | `rollerscrew`, `rollerscrew-sim-demo` |
| [`linear/`](linear/) | Piezo-hydraulic linear actuator sizing; an electroosmotic route | `linear` |
| [`zipchain/`](zipchain/) | Zip-chain push-pull actuator, limited by column buckling | `zipchain` |
| [`linear-rail-servo/`](linear-rail-servo/) | TT-gearmotor belt axis made precise by a capacitive vernier scale | `rail-servo` |
| [`corexy/`](corexy/) | The same idea in 2-D: CoreXY closed on draw-wire encoders | `corexy` |

**Motors**

| Folder | What | `make` |
|---|---|---|
| [`motor/`](motor/) | 2-D magnetostatic FEA of a 24-slot / 22-pole torque motor | `motor` |
| [`qdd-liquid/`](qdd-liquid/) | Liquid-cooled QDD motor: sizing, field solution, CAD, paper | `qdd-liquid` |
| [`pcb-motor/`](pcb-motor/) | Axial-flux PCB-stator motor, calibrated against a measured one | `pcb-motor` |
| [`chain-motor/`](chain-motor/) | Variable Chain Motor, two cells | `chain-motor` |

**Real actuators in simulation**

| Folder | What | `make` |
|---|---|---|
| [`robstride/`](robstride/) | RobStride QDDs: specs, MIT-mode model, sizing (heat decides), bench system ID, CAD | `robstride-demo` |
| [`sts3215/`](sts3215/) | SO-101 servos: stock model against fitted friction models | `sts3215-compare` |
| [`friction/`](friction/) | Extended friction models (BAM M1–M6) and a per-step MuJoCo updater | `friction-check` |
| [`mujoco/`](mujoco/) | The actuator spec bridge and test bench for the hybrid cycloidal | `actuator-bench` |

**Soft and compliant**

| Folder | What | `make` |
|---|---|---|
| [`muscle/`](muscle/) | SMA matrix muscle; HASEL cell; slow-manifold control | `muscle` |
| [`flex/`](flex/) | 2-DOF gimbal pulled by three capstan tendons | `flex` |
| [`tendril/`](tendril/) | Servo-driven TPU finger on a flexible spine | `tendril` |
| [`yzipper/`](yzipper/) | Y-zipper stiffness element (CHI '26), flexible to rigid | `yzipper` |

## How the projects are built

**One set of parameters drives everything.** The sizing model, the CAD and the simulation
read the same numbers, so they can't drift apart. Change a dimension and all three follow.

**Two layers of simulation.**
- **Analytical models answer "can this drive be efficient, strong and long-lived?"** They
  cover contact stress, friction by contact type, bearing and fatigue life, and field
  solutions for motors, and they run in milliseconds. MuJoCo can't tell you a drive's
  efficiency: give a joint a ratio and a friction loss and you've typed the efficiency in.
- **MuJoCo answers "does it work in the machine?"** Joints get the actuator's ratio,
  reflected inertia and load-dependent friction from the analytical layer, and the motors
  are current- and voltage-limited. Meshes come from the CAD.

**The CAD checks itself.** Assemblies are placed from the parameters, and each project
checks part pairs for overlap. The newer ones also check every solid is watertight and one
piece, and compare CAD mass against the sizing budget. Those checks have changed designs:
- a disc that had to be flipped rather than rotated
- bearing seats the sizing had missed
- coils that collided at their end bends

**Papers are rebuilt, then extended.** Several projects start from a paper (RV crankshaft
fatigue, the Variable Chain Motor, the Y-zipper, HASEL control). Each reproduces the paper's
own results first and states where it doesn't match.

**Renders and animation** come from the same CAD through Blender. The MuJoCo viewers run
the real kinematics.

## Setup

```bash
python3 -m venv .venv
.venv/bin/pip install build123d bd_warehouse mujoco numpy scipy matplotlib trimesh shapely
make help                      # every target
make rv                        # for example: size the RV reducer
```

Some targets need more:
- gmsh and scikit-fem for the motor field solutions (`motor/`, `qdd-liquid/`, `pcb-motor/`)
- Blender 5 for renders and animation (set `BLENDER=` if it isn't in `/opt`)
- OpenSCAD for some older previews
- Wine and FEMM 4.2 for the `qdd-liquid` cross-checks
- a display for the interactive MuJoCo viewers

Call `.venv/bin/python` directly rather than running `activate`.

## Licence

MIT, see [LICENSE](LICENSE). Vendor CAD fetched by `robstride/` isn't redistributed.
