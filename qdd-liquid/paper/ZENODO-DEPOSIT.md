# Zenodo deposit

**Concept DOI (always the latest version):** `10.5281/zenodo.23155222`

| version | DOI | file | note |
|---|---|---|---|
| 3 | `10.5281/zenodo.23166476` | `liquid-cooled-qdd-preprint-v3.pdf` | current. Field solver cross-checked against FEMM 4.2 (Table 1); jacket-only results from a thermal field solution in place of the lumped model (70% becomes 75%) |
| 2 | `10.5281/zenodo.23156072` | `liquid-cooled-qdd-preprint-v2.pdf` | superseded. Title shortened; process framing rewritten in conventional paper style (hypothesis, contributions, scope, "Design workflow" section); no change to results |
| 1 | `10.5281/zenodo.23155223` | `liquid-cooled-qdd-preprint-v1.pdf` | superseded. Record title corrected, but the PDF keeps the original title |

All published 2026-10-05. A further correction means another version on the concept DOI.
The fields used are below.

Go to <https://zenodo.org/uploads/new> and upload one file:

    qdd-liquid/paper/liquid-cooled-qdd-preprint-v3.pdf

| Field | Value |
|---|---|
| **Resource type** | Publication → **Preprint** |
| **Title** | Cooling a Quasi-Direct-Drive Joint Motor from Inside the Copper: A Simulation-Only Design Study |
| **Creator** | Newcome, Daniel · ORCID `0009-0000-6015-8713` · Affiliation: punkfab |
| **Publication date** | the day you publish |
| **License** | Creative Commons Attribution 4.0 International (CC-BY-4.0) |
| **Version** | 3 |
| **Language** | English |

**Description** (paste as-is):

> Quasi-direct-drive (QDD) robot joints trade gear ratio for transparency, and their continuous
> torque is set by how much winding heat the housing can shed. This preprint tests one
> hypothesis: that drive electronics are now efficient enough at high current that heat is the
> remaining limit, so a joint motor should be cooled directly, run hard, and geared less. The design
> studied is a Ø120 mm, 24-slot, 22-pole motor whose coils are hollow rectangular copper carrying
> dielectric oil in the bore, inside an outer jacket fed with the coldest coolant. It is evaluated
> with a first-order loss and hydraulic model, a nonlinear two-dimensional magnetostatic field
> solution that agrees with FEMM within 0.7%, a thermal field solution of the jacket-only
> alternative, and a parametric CAD assembly with interference checks. Nothing has been built or
> measured.
>
> In simulation the motor gives 16.6 N·m continuously for 850 W of heat with the copper at 126 °C,
> and 21.3 N·m at the limit of a 3 L/min pump, against about 5 N·m for the same iron cooled by air.
> The saturation limit first assumed was too pessimistic; trading copper for iron at fixed size does
> not raise torque at a given heat; and a simpler motor with solid conductors and only the jacket
> reaches about 75% of the hollow design's torque, so it is the one to build first. A 3:1 reduction
> matches a commercial 9:1 actuator's continuous torque with a ninth of the reflected inertia and
> about five times the holding power.
>
> The models and CAD were written by an AI coding agent under the author's direction; the paper
> describes that workflow and the modelling errors that cross-checks between the models exposed.
> No new cooling method is claimed.

**Keywords:** quasi-direct-drive actuator; robot joint; direct liquid cooling; hollow conductor;
permanent-magnet motor; finite element analysis; thermal management; legged robots;
AI-assisted engineering; simulation

**Related identifiers:**

| Relation | Identifier |
|---|---|
| *is supplement to* | `https://github.com/punkfab/robot-actuators` (URL) |
| *is described by* | `https://punkfab.com/blog/liquid-cooled-qdd/` (URL) |

## Before publishing

- `REFERENCES-NOTES.md` has the verification record: all 26 DOIs match Crossref, and every
  cited claim was compared with the source abstract (or the RobStride manual). Seven
  references were checked at record or title level only; the table there lists them.
- The "no published hollow-conductor robot joint motor" statement is hedged in the text:
  patents and non-English literature were not searched.
- Section 6 quotes your brief verbatim (spelling corrected) and states the AI's role.
