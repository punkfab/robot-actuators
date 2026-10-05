# Zenodo deposit

**Published 2026-10-05:** version DOI `10.5281/zenodo.23155223`, concept DOI
`10.5281/zenodo.23155222` (resolves to the latest version), <https://zenodo.org/records/23155223>.
The uploaded file is `liquid-cooled-qdd-preprint-v1.pdf` (md5 `48980eae84c4d31f1df8f3342cbfe449`).
A correction means a new version on the same concept DOI. The fields used are below.

Go to <https://zenodo.org/uploads/new> and upload one file:

    qdd-liquid/paper/liquid-cooled-qdd-preprint-v1.pdf

| Field | Value |
|---|---|
| **Resource type** | Publication → **Preprint** |
| **Title** | Cooling a Quasi-Direct-Drive Joint Motor from Inside the Copper: A Simulation-Only Design Study, Produced from a Three-Sentence Brief |
| **Creator** | Newcome, Daniel · ORCID `0009-0000-6015-8713` · Affiliation: punkfab |
| **Publication date** | the day you publish |
| **License** | Creative Commons Attribution 4.0 International (CC-BY-4.0) |
| **Version** | 1 |
| **Language** | English |

**Description** (paste as-is):

> Quasi-direct-drive (QDD) robot joints trade gear ratio for transparency, and their continuous
> torque is set by how much winding heat the housing can shed. This preprint examines one
> proposition: that drive electronics are now efficient enough at high current that heat is the
> remaining limit, so a joint motor should be cooled directly, run hard, and geared less. The design
> studied is a Ø120 mm, 24-slot, 22-pole motor whose coils are hollow rectangular copper carrying
> dielectric oil in the bore, inside an outer jacket fed with the coldest coolant. It is evaluated
> with a first-order loss and hydraulic model, a nonlinear two-dimensional magnetostatic field
> solution, a lumped thermal model of the jacket-only alternative, and a parametric CAD assembly
> with interference checks. Nothing has been built or measured.
>
> In simulation the motor gives 16.6 N·m continuously for 850 W of heat with the copper at 126 °C,
> and 21.3 N·m at the limit of a 3 L/min pump, against about 5 N·m for the same iron cooled by air.
> The saturation limit first assumed was too pessimistic; trading copper for iron at fixed size does
> not raise torque at a given heat; and a simpler motor with solid conductors and only the jacket
> reaches about 70% of the hollow design's torque, so it is the one to build first. A 3:1 reduction
> matches a commercial 9:1 actuator's continuous torque with a ninth of the reflected inertia and
> about five times the holding power.
>
> The paper also records how the study was produced: from a three-sentence brief to a paper through
> seven prompts and about seventy minutes of AI-agent working time, with the errors the models and
> checks caught in each other. No new cooling method is claimed.

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
