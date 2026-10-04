# Liquid-cooled high-current QDD motor

**Thesis under test:** solid-state drives now handle very high current efficiently, so the
limit on a joint motor is heat. Remove the heat actively, run the copper hard, and use less
gear reduction.

| file | what |
|---|---|
| `design.py` | first-order sizing: copper loss, bore hydraulics, coolant, drive (`make qdd-liquid`) |
| `fea.py` | nonlinear 2-D field solution of the same design: torque against current |
| `study.py` | the two together across geometry variants (`make qdd-liquid-fea`, about 5 minutes) |

Nothing here has been built, and there is no CAD yet.

## The concept

- **Hollow copper conductors.** Coolant flows through the bore of every turn, so the heat
  leaves from inside the copper instead of crossing insulation, iron and a housing.
- **Few turns of large conductor.** A low-voltage, high-current machine.
- **The drive on the motor,** on the same loop.
- **An outer jacket fed with the coldest coolant.** The skin stays near the loop's return
  temperature whatever the copper is doing.
- **Flow order** (regenerative in the rocket-nozzle sense: cold fluid lines the outside
  before it reaches the heat source):

      radiator -> drive cold plate -> outer jacket -> conductor bores -> radiator

- **Coolant:** dielectric oil, so it can touch live copper. Coils are electrically in series
  and hydraulically in parallel with no insulating breaks.

## Design point

Ø120 mm (RobStride RS04's diameter), 25 mm stack, 61 mm long with end turns, 24 slots / 22
poles, stator outside and bonded to the jacket, Ø54 mm free bore for a gear stage.

| | |
|---|---|
| conductor | 8 turns per coil, about 1.9 mm square, Ø1.07 mm bore |
| rated current | 70 A rms per phase at 26 A/mm² |
| continuous torque | 14.8 N·m at the motor (field solution) |
| heat at that torque | 674 W: copper 650, drive and leads 17, pump 5 |
| coolant | oil in at 60 °C, out at 73 °C, 1.85 L/min, 0.68 bar across the coils |
| temperatures | skin 60 °C, copper 93 °C at its hottest, magnets about 74 °C |
| motor mass | about 2 kg (estimate) |

## What the field solution found (`study.py`)

`design.py` on its own had to assume where the iron saturates (an air-gap shear of 90 kPa).
`fea.py` solves the field with a real B-H curve instead. The mesh is converged: torque moves
under 1% from 28k to 150k triangles.

**1. Saturation is much gentler than the assumption.**

| current density | torque | torque per amp | assumed knee said | heat |
|---|---|---|---|---|
| 14 A/mm² | 8.5 N·m | 98% | 6.8 | 188 W |
| 26 A/mm² | 14.8 N·m | 92% | 10.7 | 674 W |
| 33 A/mm² | 17.5 N·m | 86% | 12.2 | 1135 W |
| 42 A/mm² | 20.2 N·m | 77% | 13.5 | 2.0 kW |
| 100 A/mm² | 26.4 N·m | 43% | 16.6 | beyond the pump |

The air-gap field is also higher than assumed (1.04 T against 0.83 T), so every torque
figure from `design.py` alone was about 20% low before saturation and more after it.

**2. Heat and saturation arrive together.** Torque per amp is down 20% at 39 A/mm²
(19.3 N·m). A joint-sized pump (3 L/min of oil) runs out at 49 A/mm² (21.4 N·m, about
3 kW). So heat is still a limit, and the iron is the other one, at nearly the same current.

**3. More iron in the same size does not help.** Same Ø120 × 25 mm, copper traded for iron.
Torque in N·m at a given total heat:

| variant | 300 W | 850 W | 1500 W | at the cooling limit | torque per amp down 20% at |
|---|---|---|---|---|---|
| baseline: slots 50% of the pitch, 14 mm deep | 10.5 | 15.9 | 18.8 | 21.4 | 39 A/mm² |
| slots 40% | 9.7 | 15.4 | 18.8 | 22.3 | 57 A/mm² |
| slots 30% | 8.5 | 13.7 | 17.0 | 20.6 | 79 A/mm² |
| slots 10 mm deep (bigger rotor) | 10.5 | 16.7 | 20.7 | 25.0 | 65 A/mm² |
| slots 40%, 10 mm deep | 9.6 | 15.5 | 19.6 | 24.5 | 98 A/mm² |
| slots 30%, 7 mm deep | 7.7 | 12.6 | 16.1 | 18.8 | over 140 A/mm² |
| magnets 5 mm | 10.7 | 16.2 | 19.2 | 22.0 | 40 A/mm² |
| cobalt-iron laminations | 10.8 | 16.8 | 20.2 | 23.5 | 46 A/mm² |
| outrunner | 11.8 | 17.6 | 20.9 | 24.0 | 36 A/mm² |

Wider teeth do push saturation out, a long way. But they take copper away, so the same
torque needs more current density and makes more heat. The two effects cancel: every
variant lands within about 10% of the others at a given heat. The one clear gain is
shallower slots, which let the rotor grow: +5% at 850 W and +17% at the cooling limit, with
a larger free bore (Ø62 mm).

**So the trade is not size for iron; it is size for torque.** At this size, the motor makes
about 16 N·m for 850 W whatever the split between copper and iron. More torque at the same
heat means more stack length (linear) or more diameter (squared).

**4. An outrunner is worth about 11%,** not the 1.5× the first-order model suggested.

**5. What it buys at the joint.** Against RS04 (35 N·m continuous, 120 N·m peak, 9:1, about
88 W to hold rated):

| ratio | continuous output at 850 W | heat to hold 35 N·m | reflected inertia |
|---|---|---|---|
| 9:1 | 138 N·m | 43 W | 1× |
| 6:1 | 92 N·m | 95 W | 0.44× |
| 4:1 | 61 N·m | 218 W | 0.20× |
| 3:1 | 46 N·m | 415 W | 0.11× |
| 2:1 | 31 N·m | 1329 W | 0.05× |

At 3:1 the joint beats RS04's continuous torque with a ninth of the reflected inertia, and
spends about five times the power holding 35 N·m. Direct drive is out of reach at this
size: 20 N·m at the cooling limit.

![field-solved study](out/study.png)

## What the sizing model found (`design.py`)

These hold with the field solution's torque numbers.

**"Very high current" is a choice of turns, not a goal.** What matters is current density.
Fewer turns raise the amps and the drive loss without adding torque: at the same torque,
2 turns per coil means about 280 A and 195 W in the switches and leads, 8 turns means 70 A
and 17 W. More turns shrink the bore, which costs pump pressure. `study.py` picks 8.

**Running the loop hot costs torque and saves radiator.** Copper resistance rises 0.39% per
kelvin and the magnets weaken. Going from a 60 °C to a 100 °C coolant return adds about 12%
heat at the same current and needs a higher magnet grade; the payoff is a radiator a bit
over half the size.

**The cold-first flow order works.** The skin sits within a degree or two of the coolant
return temperature at every operating point.

![sizing model](out/design.png)

## Where this leaves the thesis

Mostly intact. Cooling the copper from inside lifts continuous torque from about 5 N·m
(the same iron, air-cooled, 88 W) to 15–21 N·m, and lets the gear ratio drop from 9:1 to
about 3:1. Two corrections:

- Heat is not the only limit. The iron runs out at nearly the same current the pump does.
- The price is energy. Torque rises with current and heat with its square, so a low-ratio
  joint holding a load burns several times what a geared one does.

## Not modelled

- 3-D effects: the field solution is 2-D, and a 25 mm stack will lose some torque to end
  leakage.
- Torque ripple (one rotor position), demagnetization, iron loss beyond a rough estimate.
- The B-H curve is a generic electrical-steel fit; cobalt-iron is that curve stretched 15%.
- AC loss, beyond noting the skin depth (about 10 mm) is far larger than the conductor.
- Coil bends, manifolds and fittings in the pressure drop; flow sharing between coils.
- Making it: a Ø1 mm bore in a 1.9 mm conductor, its joints, seals and insulation.
- Coolant properties are handbook approximations. Magnet temperature is an estimate.
