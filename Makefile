# Robot actuators — build (CAD -> STEP/STL) and sim (kinematic MuJoCo viewer) targets.
# All parts use the in-repo virtualenv (build123d + mujoco live there).
PY := .venv/bin/python

.DEFAULT_GOAL := help

.PHONY: help
help:
	@echo "targets:"
	@echo "  make cycloidal-center       build the straddle-carrier centre-output cycloidal"
	@echo "  make sim-cycloidal-center   build + open its kinematic MuJoCo viewer"
	@echo "  make planetary              build the standalone planetary reducer"
	@echo "  make sim-planetary          build + open its kinematic MuJoCo viewer"
	@echo "  make planetary-inverted     build the inverted carrier-out planetary (NEMA-17 in/out)"
	@echo "  make sim-planetary-inverted build + open its kinematic MuJoCo viewer"
	@echo "  make cycloidal-2stage       build the two-stage compound cycloidal (body output, yoke)"
	@echo "  make sim-cycloidal-2stage   build + open its kinematic MuJoCo viewer (121:1)"
	@echo "  make cycloidal              build the compound planetary+cycloidal drive"
	@echo "  make motor                  EM FEA: mesh + flux-density field of the torque motor"
	@echo "  make motor-curves           EM FEA: + cogging & Kt sweeps (slower)"
	@echo "  make linear                 piezo-hydraulic linear actuator sizing + F-v envelope"
	@echo "  make linear-eo              add the electroosmotic route + piezo-vs-EO comparison"
	@echo "  make rail-encoder           capacitive vernier scale (digital-caliper) resolution model"
	@echo "  make rail-servo             control proof: TT gearmotor + cap scale cancels backlash"
	@echo "  make rail-cad               build the linear rail servo CAD (STEP/STL parts)"
	@echo "  make sim-rail               build CAD + open its kinematic MuJoCo viewer"
	@echo "  make drawwire               draw-wire (tape-measure) encoder resolution model"
	@echo "  make corexy                 CoreXY proof: draw-wire loop cancels belt backlash/stretch"
	@echo "  make pcb-motor              axial-flux PCB-stator motor model + reduction-need calculator"
	@echo "  make pcb-motor-fea          2D unrolled magnetostatic FEA cross-check of the air-gap flux"
	@echo "  make pcb-motor-benchmark    calibrate the motor model vs a measured PCB motor (Wang 2025)"
	@echo "  make flex                   2-DOF tendon-gimbal statics: 3 capstans, tensions, workspace"
	@echo "  make flex-cad               build the tendon-gimbal CAD (STEP/STL parts)"
	@echo "  make sim-flex               build CAD + open its 2-DOF MuJoCo viewer"
	@echo "  make tendril                servo-driven TPU continuum finger: curl/force statics"
	@echo "  make tendril-cad            build the tendril CAD (TPU finger + PLA servo mount)"
	@echo "  make sim-tendril            build CAD + open its curling MuJoCo viewer"
	@echo "  make zipchain               zip-chain LINEAR actuator: buckling-limited push envelope"
	@echo "  make zipchain-cad           build the zip-chain CAD (2 strands + sprocket + head)"
	@echo "  make sim-zipchain           build CAD + open its deploy MuJoCo viewer"
	@echo "  make yzipper                Y-zipper 3-strip flex<->rigid: stiffness-switch model"
	@echo "  make yzipper-cad            build the Y-zipper CAD (3 strips + 3-way slider)"
	@echo "  make sim-yzipper            build CAD + open its soft->rigid MuJoCo viewer"
	@echo "  make chain-motor            2-cell Variable Chain Motor: kinematics, 2S/2P envelopes, eff."
	@echo "  make chain-motor-cad        build its CAD (gears, swing link) + interference check"
	@echo "  make chain-motor-check      MuJoCo dynamic checks (hinge neutrality, lead swap, inertia)"
	@echo "  make sim-chain-motor        build CAD + open its bending-chain MuJoCo viewer"
	@echo "  make friction               extended friction models (BAM M1-M6): fitted-servo η table"
	@echo "  make friction-check         MuJoCo pendulum bench: sim obeys the friction models"
	@echo "  make actuator-bench         hybrid 40:1 on the MuJoCo test bench (scenarios A-F)"
	@echo "  make sts3215-compare        SO-101: stock MJCF servo vs BAM STS3215 models (sag, step...)"
	@echo "  make sts3215-crosscheck     re-run software-mfg's workcell gate with BAM servos"
	@echo "  make sts3215-synthetic      replay + fit pipeline self-test on synthetic logs"
	@echo "  make robstride              RobStride spec table + thermal fit to the overload tables"
	@echo "  make robstride-demo         size, tune and verify a 3-DOF arm on RobStrides"
	@echo "  make robstride-bench        pendulum-bench sysid self-test (recover a known RS03)"
	@echo "  make robstride-cad          fetch vendor STEPs, normalize (output flange up), line-up"
	@echo "  make qdd-liquid             direct-liquid-cooled high-current QDD motor: first-order sizing"
	@echo "  make qdd-liquid-fea         nonlinear field solution + copper/iron variants (about 5 minutes)"
	@echo "  make qdd-liquid-cad         shallow-slot motor assembly: STEP, interference checks, cutaway render"
	@echo "  make qdd-liquid-jacket      same motor, solid conductors, jacket-only cooling: what it gives up"
	@echo "  make qdd-liquid-render      Blender renders of the motor: cutaway and wound stator"
	@echo "  make qdd-liquid-femm        cross-check of the field solver against FEMM 4.2 (needs wine + FEMM)"
	@echo "  make qdd-liquid-thermal     jacket-only motor: thermal field solution in FEMM against the lumped model"
	@echo "  make qdd-liquid-paper       figures and numbers for the paper (qdd-liquid/paper); -paper-pdf builds the PDF"
	@echo "  make rollerscrew            roller-screw linear actuator (8 kN knee class): sizing against a day's duty"
	@echo "  make rollerscrew-cad       its CAD assembly: STEP, STLs, interference checks"
	@echo "  make rollerscrew-render    Blender renders: half-section and the roller screw alone"
	@echo "  make rollerscrew-sim       MuJoCo leg with the actuator in its knee: interactive viewer (sliders)"
	@echo "  make rollerscrew-sim-demo  statics check, landing A/B across leads, loaded squat -> out/sim.png, landing.gif"
	@echo "  make rollerscrew-sim-selftest  build + settle + one landing, no window"
	@echo "  make rv-rcf                 RV crankshaft rolling-contact fatigue (Li et al. 2022) against the paper"
	@echo "  make rv                     RV reducer actuator: kinematics, loads, crank life, motor"
	@echo "  make rv-cad                 RV actuator CAD: STEP, interference and part checks (~9 min)"
	@echo "  make rv-render              Blender half-section through a crankshaft"
	@echo "  make rv-sim                 MuJoCo kinematic playback of the RV actuator (interactive)"
	@echo "  make rv-printed             3D-printable RV: ball bearings + steel dowels, load ranking, CAD, BOM"
	@echo "  make rv-explode             exploded drawing of the printable RV with balloons and a parts list"
	@echo "  make rv-assemble-video      Blender animation: parts fly together, then the reducer runs (~45 min; --quick)"
	@echo "  make rv-sim-video           meshing check + render out/rv_sim.mp4 / .gif"
	@echo "  make muscle                 matrix muscle: SMA cell, bundle, crosstalk sheet"
	@echo "  make muscle-sm              slow-manifold feed-forward (Bettini 2026) on the SMA cell"
	@echo "  make muscle-hasel           HASEL pouch cell (Kellaris model), sized beside the SMA unit"
	@echo "  make muscle-compare         SMA vs HASEL on the same control test (a few minutes)"

# --- centre-output cycloidal ------------------------------------------------
.PHONY: cycloidal-center
cycloidal-center:
	$(PY) cycloidal-center/drive.py

.PHONY: sim-cycloidal-center
sim-cycloidal-center: cycloidal-center
	$(PY) cycloidal-center/sim.py $(REV)

# --- standalone planetary reducer -------------------------------------------
.PHONY: planetary
planetary:
	$(PY) planetary/reducer.py

.PHONY: sim-planetary
sim-planetary: planetary
	$(PY) planetary/sim.py $(REV)

# --- inverted (carrier-out) planetary, NEMA-17 in/out -----------------------
.PHONY: planetary-inverted
planetary-inverted:
	$(PY) planetary-inverted/drive.py

.PHONY: sim-planetary-inverted
sim-planetary-inverted: planetary-inverted
	$(PY) planetary-inverted/sim.py $(REV)

# --- compound planetary+cycloidal drive -------------------------------------
.PHONY: cycloidal
cycloidal:
	$(PY) cycloidal/drive.py

# --- two-stage compound cycloidal (body output, external-yoke straddle) ------
.PHONY: cycloidal-2stage
cycloidal-2stage:
	$(PY) cycloidal-2stage/drive.py

.PHONY: sim-cycloidal-2stage
sim-cycloidal-2stage: cycloidal-2stage
	$(PY) cycloidal-2stage/sim.py $(REV)

# --- custom torque-motor electromagnetic FEA (gmsh + scikit-fem) ------------
.PHONY: motor
motor:
	$(PY) motor/fea.py

.PHONY: motor-curves
motor-curves:
	$(PY) motor/fea.py curves

# --- piezo-hydraulic linear actuator (first-order sizing model) -------------
.PHONY: linear
linear:
	$(PY) linear/actuator.py

.PHONY: linear-eo
linear-eo:
	$(PY) linear/actuator.py eo

# --- linear rail servo: crappy TT gearmotor + capacitive vernier scale -------
.PHONY: rail-encoder
rail-encoder:
	$(PY) linear-rail-servo/encoder.py

.PHONY: rail-servo
rail-servo:
	$(PY) linear-rail-servo/servo.py

.PHONY: rail-cad
rail-cad:
	$(PY) linear-rail-servo/rail.py

.PHONY: sim-rail
sim-rail: rail-cad
	$(PY) linear-rail-servo/sim.py $(REV)

# --- CoreXY stage on TT motors, closed on draw-wire (string-pot) encoders -----
.PHONY: drawwire
drawwire:
	$(PY) corexy/drawwire.py

.PHONY: corexy
corexy:
	$(PY) corexy/corexy.py

# --- axial-flux PCB-stator motor (analytical model + reduction calculator) ----
.PHONY: pcb-motor
pcb-motor:
	$(PY) pcb-motor/motor.py

.PHONY: pcb-motor-fea
pcb-motor-fea:
	$(PY) pcb-motor/fea.py

.PHONY: pcb-motor-benchmark
pcb-motor-benchmark:
	$(PY) pcb-motor/benchmark.py

# --- 2-DOF tendon-driven gimbal (3 capstans, sprung center) ------------------
.PHONY: flex
flex:
	$(PY) flex/flex.py

.PHONY: flex-cad
flex-cad:
	$(PY) flex/cad.py

.PHONY: sim-flex
sim-flex: flex-cad
	$(PY) flex/sim.py $(REV)

# --- servo-driven continuum tendril (1 hobby servo, antagonistic strings) ------
.PHONY: tendril
tendril:
	$(PY) tendril/tendril.py

.PHONY: tendril-cad
tendril-cad:
	$(PY) tendril/cad.py

.PHONY: sim-tendril
sim-tendril: tendril-cad
	$(PY) tendril/sim.py $(REV)

# --- zip-chain LINEAR actuator (2-strip rigid chain, buckling-limited push) -----
.PHONY: zipchain
zipchain:
	$(PY) zipchain/zipchain.py

.PHONY: zipchain-cad
zipchain-cad:
	$(PY) zipchain/cad.py

.PHONY: sim-zipchain
sim-zipchain: zipchain-cad
	$(PY) zipchain/sim.py $(REV)

# --- Y-zipper 3-strip flexible<->rigid STIFFNESS element (closed Δ section) ------
.PHONY: yzipper
yzipper:
	$(PY) yzipper/yzipper.py

.PHONY: yzipper-cad
yzipper-cad:
	$(PY) yzipper/cad.py

.PHONY: sim-yzipper
sim-yzipper: yzipper-cad
	$(PY) yzipper/sim.py $(REV)

# --- 2-cell Variable Chain Motor (Tada et al. IROS 2025): geared BLDC cells that bend ---
.PHONY: chain-motor
chain-motor:
	$(PY) chain-motor/chain.py

.PHONY: chain-motor-cad
chain-motor-cad:
	$(PY) chain-motor/cad.py

.PHONY: chain-motor-check
chain-motor-check: chain-motor-cad
	$(PY) chain-motor/sim.py check

.PHONY: sim-chain-motor
sim-chain-motor: chain-motor-cad
	$(PY) chain-motor/sim.py $(REV)

# --- extended gearbox friction (after Rhoban BAM, ICRA 2025) -----------------
.PHONY: friction
friction:
	$(PY) friction/models.py

.PHONY: friction-check
friction-check:
	$(PY) friction/bench.py

.PHONY: actuator-bench
actuator-bench:
	$(PY) mujoco/run.py

# --- STS3215 / SO-101 real2sim (BAM actuator model) --------------------------
.PHONY: sts3215-compare sts3215-crosscheck sts3215-synthetic
sts3215-compare:
	$(PY) sts3215/compare.py

sts3215-crosscheck:
	$(PY) sts3215/crosscheck.py

sts3215-synthetic:
	$(PY) sts3215/replay.py --synthetic
	$(PY) sts3215/fit.py --synthetic

# --- RobStride QDD actuators ------------------------------------------------
.PHONY: robstride robstride-demo robstride-bench robstride-cad
robstride:
	$(PY) robstride/specs.py
	$(PY) robstride/thermal.py
	$(PY) robstride/can_mit.py

robstride-demo:
	$(PY) robstride/demo.py

robstride-bench:
	$(PY) robstride/bench.py --synthetic RS03

robstride-cad:
	$(PY) robstride/cad.py

.PHONY: muscle
muscle:
	$(PY) muscle/cell.py
	cd muscle && ../$(PY) bundle.py && ../$(PY) sheet.py

.PHONY: muscle-sm
muscle-sm:
	$(PY) muscle/slowmanifold.py

.PHONY: muscle-hasel
muscle-hasel:
	$(PY) muscle/hasel.py

.PHONY: muscle-compare
muscle-compare:
	$(PY) muscle/compare.py

.PHONY: qdd-liquid
qdd-liquid:
	$(PY) qdd-liquid/design.py

.PHONY: qdd-liquid-fea
qdd-liquid-fea:
	$(PY) qdd-liquid/fea.py
	$(PY) qdd-liquid/study.py

.PHONY: qdd-liquid-cad
qdd-liquid-cad:
	$(PY) qdd-liquid/cad.py

.PHONY: qdd-liquid-jacket
qdd-liquid-jacket:
	$(PY) qdd-liquid/jacket.py

.PHONY: qdd-liquid-render
qdd-liquid-render:
	$(PY) qdd-liquid/render.py

# FEMM 4.2 under wine (WINEPREFIX=~/.wine-femm, C:\\femm42); see qdd-liquid/README.md
.PHONY: qdd-liquid-femm
qdd-liquid-femm:
	$(PY) qdd-liquid/femm_check.py

.PHONY: qdd-liquid-thermal
qdd-liquid-thermal:
	$(PY) qdd-liquid/thermal_check.py

.PHONY: qdd-liquid-paper
qdd-liquid-paper:
	$(PY) qdd-liquid/paper/figures.py

# tectonic fetches the LaTeX packages it needs; override with TECTONIC=/path/to/tectonic
TECTONIC ?= tectonic
.PHONY: qdd-liquid-paper-pdf
qdd-liquid-paper-pdf:
	cd qdd-liquid/paper && $(TECTONIC) main.tex
	@echo "wrote qdd-liquid/paper/main.pdf"

# --- roller-screw linear actuator (Optimus-style) ----------------------------
.PHONY: rollerscrew
rollerscrew:
	$(PY) rollerscrew/design.py

.PHONY: rollerscrew-cad
rollerscrew-cad: rollerscrew
	$(PY) rollerscrew/cad.py

.PHONY: rollerscrew-render
rollerscrew-render: rollerscrew-cad
	$(PY) rollerscrew/render.py

.PHONY: rollerscrew-sim
rollerscrew-sim:
	$(PY) rollerscrew/sim.py

.PHONY: rollerscrew-sim-demo
rollerscrew-sim-demo:
	$(PY) rollerscrew/sim.py --demo

.PHONY: rollerscrew-sim-selftest
rollerscrew-sim-selftest:
	$(PY) rollerscrew/sim.py --selftest

# --- RV reducer actuator (after Li et al. 2022) -------------------------------
.PHONY: rv-rcf
rv-rcf:
	$(PY) rv/rcf.py

.PHONY: rv
rv:
	$(PY) rv/design.py

.PHONY: rv-cad
rv-cad: rv
	$(PY) rv/cad.py

.PHONY: rv-render
rv-render: rv-cad
	$(PY) rv/render.py

.PHONY: rv-sim
rv-sim:
	$(PY) rv/sim.py

.PHONY: rv-sim-video
rv-sim-video:
	$(PY) rv/sim.py --video

.PHONY: rv-printed
rv-printed:
	$(PY) rv/printed.py

.PHONY: rv-explode
rv-explode:
	$(PY) rv/explode.py

.PHONY: rv-assemble-video
rv-assemble-video:
	$(PY) rv/assemble_anim.py
