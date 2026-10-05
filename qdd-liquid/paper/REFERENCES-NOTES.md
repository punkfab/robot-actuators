# Reference verification notes (qdd-liquid preprint)

Verified 2026-10-04. BibTeX is in `references.bib`.

## Summary

All 15 candidates were confirmed to exist, except that item 15 needed a substitute source. No hollow-conductor, direct-cooled motor for a robot joint turned up, but direct winding cooling in a legged-robot actuator (immersion) is already published, so the novelty claim has to be narrow.

Bibliographic data comes from Crossref records, plus MIT DSpace (Katz) and the Aaltodoc repository record (Arkkio). Relevance notes come from publisher abstracts via Crossref/OpenAlex unless marked otherwise. No full texts were opened, and IEEE Xplore and Wiley pages blocked fetching.

## Part 1: verification notes

1. **wensing2017: confirmed.** Introduces the "impact mitigation factor" for backdrivability; the MIT Cheetah leg controls contact forces in bounding with contact times down to 85 ms and peak forces over 450 N.
2. **seok2015: confirmed.** IEEE/ASME TMECH 20(3):1117-1129, June 2015 issue (online 2014). The 33 kg robot trots at 6 m/s on 973 W, cost of transport 0.5; 76% of total energy is motor heat loss, which is a useful motivating number.
3. **katz2018: confirmed.** Benjamin G. Katz, S.M. thesis, MIT Dept. of Mechanical Engineering, 2018, advisor Sangbae Kim. Low-cost modular actuator and a 12-DOF quadruped that does a 360 degree backflip; the abstract gives no thermal numbers.
4. **urata2010: confirmed.** IROS 2010, pp. 4497-4502. The abstract describes "active temperature control" combining motor internal-temperature estimation with forced liquid cooling, plus high-power drivers; it gives no numbers. The often-quoted "over 1000 deg/s, 350 N·m knee" figures appeared only in secondary sources, so check the paper before quoting them. A same-titled Japanese journal version also exists (J. Robotics Society of Japan 28(7):865-871, 2010, doi:10.7210/jrsj.28.865).
5. **paine2015: confirmed.** Actuators 4(3):182-202. Retrofitted liquid cooling of an off-the-shelf motor (housing cooling, not winding-direct): continuous torque x2.58 (within 9% of datasheet expectation), continuous power x2, with 2.2% lower efficiency than air cooling.
6. **kozuki2016: confirmed.** IROS 2016, pp. 2135-2140. Only the first abstract sentence could be retrieved: the skeletal structure is used to release heat by latent heat as well as carry load. No numbers verified.
   - Author-name caution: Crossref lists the second author as "Hirose Toshinori" (Semantic Scholar: Toshinori Hirose) and spells "Shinske Nakashima" (elsewhere "Shinsuke"). The BibTeX uses the normalised spellings.
7. **Lindh: both papers confirmed.**
   - **lindh2016** (IEEE TEC 31(4):1257-1266): 100 kW axial-flux double-stator machine with helical tooth coils of a hybrid conductor (stainless-steel coolant tube wrapped in Litz wire), polyalphaolefin oil coolant; reported feasible with significant thermal improvement.
   - **lindh2017** (IEEE TIE 64(8):6086-6095): the same hybrid conductor on a tooth-coil axial-flux machine, compared against water-jacket indirect cooling, with oil or water; practical "also in small machines", best where stator Joule losses dominate.
   - Correction: the coolant runs in a stainless tube inside the Litz bundle, not in a hollow copper conductor. The abstracts give no gain factors.
8. **schiefer2015: confirmed.** IEMDC 2015, pp. 1820-1825. Water channels in the unused slot space of flat-wire concentrated windings, compared by FEM/CFD and single-tooth measurements against a water-jacket round-wire IPMSM; the jacket can be omitted. No numeric gain in the abstract.
9. **Hollow conductors: three confirmed.**
   - **wohlers2018** (Electrical Engineering 100(4):2299-2308): direct-liquid-cooled tooth-coil designs based on cast coils with slot fill up to 90%, possible current density 100 A/mm², additional losses cut to about 50%. These numbers are from an abstract snippet in search results; Crossref had no abstract.
   - **wu2020** (ICEM 2020, pp. 1497-1503): additively manufactured (DMLS) aluminium-alloy hollow conductors with coolant flow in a direct-drive concentrated-winding PM motor; targets above 20 kW/kg and above 95% efficiency.
   - **wu2021** (IEEE TIA 57(4):3632-3642): hollow conductors integrated with heat pipes, 250 kW aircraft machine, AlSi10Mg coil samples printed.
10. **Reviews: both confirmed.**
    - **gai2019** (IEEE TIE 66(3):1681-1692): review of air, liquid and phase-change cooling schemes with real traction-motor examples and computation methods.
    - **yang2017** (IET Electrical Systems in Transportation 7(2):104-116, online 2016): review of losses, materials, cooling techniques and thermal analysis.
    - Neither abstract names direct conductor cooling explicitly; check the body before citing for that specific point.
11. **brauer1975: confirmed.** IEEE Trans. Magn. MAG-11(1):81, a single page. Fits B-H curves as H = (k1·exp(k2·B²) + k3)·B, with simple reluctivity and derivative expressions for finite elements.
12. **arkkio1987: confirmed.** Acta Polytechnica Scandinavica, Electrical Engineering Series No. 59, Helsinki University of Technology, 1987, 97 pp., ISBN 951-22-6076-X, ISSN 0001-6845. Combined 2D finite-element field and winding circuit equations for induction motors.
    - The title has no "the" before "induction motors" (some citing papers add one).
    - A search snippet gave a different ISBN (951-666-250-1), probably the print series edition; not verified.
13. **geuzaine2009: confirmed.** IJNME 79(11):1309-1331, 2009, doi:10.1002/nme.2579.
14. **pyrhonen2013: confirmed.** 2nd edition, Wiley, 2013, ISBN 978-1-118-58157-5, doi:10.1002/9781118701591. The 1st edition is 2008, doi:10.1002/9780470740095.
15. **Turbo-generator precedent: partly confirmed.**
    - **klempner2018** (Klempner and Kerszenbaum, Handbook of Large Turbo-Generator Operation and Maintenance, 3rd ed., Wiley-IEEE Press, 2018; Ch. 2 "Generator Design and Construction", pp. 53-168) exists and its chapter abstract discusses directly cooled stator windings. The hollow-strand sentence itself was not read (Wiley returned 403), so check the page before citing.
    - **bauer2021** (ASME Power 2021) has a citable sentence that was read: direct cooling is needed above about 250 MVA, with the medium "introduced via hollow conductors within the stator bars"; above about 700 MVA direct water cooling is used; the abstract refers to 50 years of practice. Crossref gives no page or paper number.
    - **svoboda2004** (ASME Power 2004, pp. 373-377): "Copper is a traditional material for hollow conductors used in generator stator water cooling", with 30 years of experience on stainless hollow conductors.

## Part 2: robot-specific prior art

**No published hollow-conductor (coolant inside the winding conductor) motor for a robot joint, legged, humanoid or QDD, was found.** However, direct winding cooling of a legged-robot actuator already exists (hit 1), and in-slot cooling of a humanoid joint motor exists (hit 2). A defensible claim is "hollow-conductor direct cooling applied to a QDD joint motor (simulation study)", not "first direct-cooled robot actuator".

1. **zhu2021** (Zhu, Ahn, Hong, ICAR 2021, pp. 1137-1143): integrated direct coil immersion cooling of a BLDC motor for legged robots, built and tested. Thermal resistance 3.7x better than liquid (housing) cooling and 10x better than forced air; in some cases the continuous-to-peak torque gap is "effectively eliminated"; also compares two-phase with single-phase cooling. This is the closest prior art.
2. **li2026** (Li, Zhou, Mei, Wang, Applied Thermal Engineering 298, 131111, 2026): hydrogel phase-change cooling structures embedded in the stator slots of a humanoid joint motor. Peak temperature 24.5 °C lower and operating time 43.1% longer; these numbers are from a publisher-abstract snippet in search results.
3. **zhu2019** (Zhu, Hooks, Hong, IEEE/ASME AIM 2019, pp. 36-43): liquid-cooled proprioceptive actuator (coolant channels in the housing) with a multi-actuator thermal model for peak versus continuous torque trade-offs. The figures "32 N·m peak, about 8 N·m continuous on air, about 21 N·m with liquid" come from a secondary description, not the abstract.
4. **mazumdar2016** (Mazumdar et al., Sandia, ASME DSCC 2016, paper V002T26A004): finned BLDC housing gives 50% more heat transfer than nominal; 79% more with forced air and 107% more with pumped water; net energy savings of 4-6%.
5. **kim2018** (Kim, Ahn, Campbell, Paine, Sentis, IEEE/ASME TMECH 23(6):2704-2714): viscoelastic liquid-cooled actuator, with liquid cooling at the motor housing and a ball-screw series-elastic drive; the abstract has no cooling numbers.

Two further items, not counted among the five:

- **geelen2024** (IEEE Trans. Magn. 60(9), 2024): simulation study of direct-water-cooled hollow conductors in coreless linear motors. It is a precision motion system, not a robot joint, but it is a small-machine hollow-conductor modelling study worth citing.
- **jin2025** (Jin, Kobayashi, Doi, Humanoids 2025, pp. 1055-1062, doi:10.1109/humanoids65713.2025.11203109): "Development of the Leg Structure with Concentrated Liquid-Cooling Actuators for Long-Term Operation". Existence confirmed in Crossref only; the abstract could not be read, so the cooling type is unknown. It is not in the BibTeX.

**What was searched:** Crossref and OpenAlex queries, plus web searches, combining hollow conductor / direct winding cooling / oil immersion / spray / in-slot cooling with robot joint / legged / humanoid / QDD / exoskeleton / proprioceptive actuator.

**What was not searched:** Google Scholar directly, patents, or Chinese/Japanese/Korean-language literature. A patent or non-English hollow-conductor robot motor cannot be ruled out.

## Notes on the BibTeX

- kozuki2016 uses the normalised author spellings "Hirose, Toshinori" and "Nakashima, Shinsuke"; the Crossref/IEEE record differs (see item 6).
- mazumdar2016's booktitle is a shortened form of the long ASME volume title; the paper number, DOI and authors are from Crossref.
- klempner2018's address "Hoboken, NJ" and pyrhonen2013's address "Chichester, UK" are standard Wiley imprint locations that were not verified against the title pages.

## Not confirmed or corrected

- **Nothing in items 1-14 failed to verify.**
- **Item 7 correction:** the Lindh conductor is a stainless coolant tube wrapped in Litz wire, not a hollow copper conductor.
- **Item 12 correction:** exact title has no "the" before "induction motors".
- **Item 15, klempner2018:** book and chapter exist, but the hollow-strand sentence was not read; bauer2021 and svoboda2004 have abstract text that was read.
- **Numbers not from primary abstracts (check before quoting):** Urata "1000 deg/s, 350 N·m"; Zhu 2019 "8 to 21 N·m continuous, 32 N·m peak"; Wohlers "100 A/mm², 90% fill, about 50% additional losses"; Li 2026 "24.5 °C, 43.1%".
- **Abstracts not obtained in full:** kozuki2016 (first sentence only) and jin2025 (none).
- **gai2019 and yang2017:** whether they cover direct conductor cooling specifically is not confirmed from the abstracts.


## Second pass, 2026-10-05: every claim in the paper against its source

Two checks were run on the paper as built.

**1. Records.** Every DOI in `references.bib` was fetched from the Crossref API and its
title, authors, year, volume, issue and pages compared with the entry. All 26 match. The
two theses were checked against their repository records (MIT DSpace, Aaltodoc), the two
textbooks against Open Library by ISBN, and the two repository URLs resolve.

**2. Claims.** Each sentence that cites a source was compared with that source's abstract
(Crossref or OpenAlex), or with the document itself for the RobStride manual. No full
texts of paywalled papers were opened.

| key | what the paper says | checked against | result |
|---|---|---|---|
| seok2015, wensing2017 | QDD: high-torque-density motor, single low-ratio stage, backdrivable, force control at the joint | abstracts | supported |
| seok2015 | 76% of energy in trotting is motor heat | abstract | supported, exact figure |
| katz2018 | low-cost modular actuator | DSpace abstract | supported |
| urata2010 | forced liquid cooling with motor internal-temperature estimation | abstract | supported; wording corrected from "winding-temperature" |
| paine2015 | retrofitted liquid cooling, 2.58x continuous torque, measured | abstract | supported, exact figure |
| kim2018, zhu2019 | built and modelled liquid-cooled actuators | abstracts | supported |
| mazumdar2016 | redesigned housing +50% heat transfer, +107% with pumped water | abstract | supported; "finned" removed (not in the abstract) |
| kozuki2016 | skeleton releases heat by latent heat | title and first abstract sentence | supported at that level only |
| li2026 | evaporative cooling of humanoid joint motors with phase-change hydrogels | title only (no abstract available) | reduced to what the title says; "in the stator slots" removed |
| zhu2021 | coil immersion, built and tested, 3.7x vs liquid cooling, 10x vs forced air | abstract | supported, exact figures; "housing" removed |
| bauer2021, svoboda2004 | hollow water-cooled stator conductors are long-standing in large turbo-generators | abstracts | supported |
| lindh2016, lindh2017 | axial-flux machine, stainless tube wrapped in Litz wire, oil or water, compared with water jackets | abstracts | supported; "oil" corrected to "oil or water" |
| schiefer2015 | water channels in unused slot space of a concentrated winding | abstract | supported |
| wohlers2018 | direct liquid cooling of tooth-coil windings | title only (no abstract available) | reduced to what the title says; "cast" removed |
| wu2020, wu2021 | additively manufactured hollow conductors for aircraft machines | abstracts | supported |
| geelen2024 | simulation study, water-cooled hollow conductors, coreless linear motors | abstract | supported; "small" removed |
| gai2019, yang2017 | reviews of machine cooling | abstracts | supported |
| brauer1975 | reluctivity fit k1 exp(k2 B^2) + k3 | abstract | supported, same form |
| geuzaine2009 | Gmsh mesh generator | abstract | supported |
| arkkio1987 | air-gap torque integral ("Arkkio's method") | repository record only | record verified; the formula is the one conventionally attributed to this thesis, the thesis text was not opened |
| pyrhonen2013 | standard machine sizing relations | book record | record verified; no page checked |
| bergman2011 | Nu = 4.36 laminar, Dittus-Boelter turbulent | Open Library record | record verified; standard textbook content, no page checked |
| sutton2017 | regenerative cooling: propellant lines the nozzle wall before combustion | Open Library record | record verified; standard textbook content, no page checked |
| robstride-rs04 | 9:1, 35 N·m rated, 120 N·m peak, 1.42 kg, 42.2 mm, Kt 2.1 N·m/Arms, 0.16 ohm line | RS04 User Manual 260713 and spec booklet text | all found in the documents; the 35 N·m rating is on a 220 x 200 mm heat sink, now stated in the paper |

**Changes made to the paper in this pass:** six sentences reworded as noted above, the
RS04 heat-sink condition added, two textbook references added for the heat-transfer
correlations and the regenerative-cooling analogy, and the uncited `klempner2018` entry
removed from the bib.

**Still not verified against full text:** arkkio1987, pyrhonen2013, bergman2011 and
sutton2017 (records only), and kozuki2016, li2026 and wohlers2018 (title-level claims
only). The paper's statements about these are limited accordingly.
