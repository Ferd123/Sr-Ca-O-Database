# Ca-Sr-O thermodynamic database for Thermo-Calc

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21925958.svg)](https://doi.org/10.5281/zenodo.21925958)

CALPHAD description of the Ca-Sr-O system with the liquid phase in the **associate
model** (species Ca, Sr, O, CaO, SrO), assembled from five published sources that used
mutually incompatible models.

**Author:** Fernando Ivan Cruz-Cervantes (Instituto Politécnico Nacional)
Assembled: 2026-08-03. Thermo-Calc 2025a.
License: CC BY 4.0 (database and data) / MIT (code). See section 10.

---

## 1. Files

| File | What it is |
|---|---|
| `CaSrO.TDB` | The database. Condensed phases + partial GAS. |
| `CaSrO_opt.tdb`, `CaSrO_opt_liq.tdb` | Variants with the liquid reoptimized. |
| `validar_CaSrO.TCM` | Console Mode macro, 7 blocks (see below). |
| `CaSrO.POP`, `optimizar_CaSrO_PARROT.TCM` | Optimization with PARROT. |
| `experimental_CaSrO.EXP` | Experimental data digitized from the literature, with the provenance of each dataset in the header. |
| `*.csv` | Tabulated experimental points and the calculated diagrams in CSV. |
| `dos.exp`, `diagrama tc.exp` | **Calculated** diagrams, exported from Console Mode in DATAPLOT format. These are results, not experimental data. |
| `*.py` | Assessment and figure scripts. |

The two calculated `.exp` files correspond to the two configurations of section 7:
`diagrama tc.exp` comes from the run with all 7 phases active (`LIQUID`, `HALITE#1/#2`,
`FCC_A1`, `BCC_A2`, `SRO2` appear), and `dos.exp` from the pseudobinary section with
everything suspended except liquid and halite (only `LIQUID` and `HALITE#1/#2`). The
plotting scripts read `dos.exp`.

The papers the parameters come from are **not** in the repository: they are under
copyright. Full citations are in section 2 and in the `LIST_OF_REFERENCES` of the TDB
itself.

The macro blocks, separated by `@&` so they can be run one at a time:

| # | What it does | What should come out |
|---|---|---|
| 1 | Loads the database and lists the system | 7 phases |
| 2 | Congruent melting of CaO | 3222.0 K |
| 3 | Congruent melting of SrO | 2870.5 K |
| 4 | Ca-O diagram (GAS suspended) | eutectic ~1094 K, 1.7 at.%O |
| 5 | CaO-SrO pseudobinary | gap x(Sr) 0.116–0.348 at 1100 K |
| 6 | Metallic Ca-Sr binary | continuous fcc and bcc, no intermediate phases |
| 7 | Stability of CaO₂ against pO₂ | decomposition at low T |

---

## 2. Sources

| Ref | Citation | What it contributes |
|---|---|---|
| [1] | Deffrennes, Jakse, Alvares, Nuta, Pasturel, Khvan, Pisch, *Calphad* **69** (2020) 101764 | CaO(cr), CaO₂(cr), CaO(liq), L(Ca,CaO), L(CaO,O) |
| [2] | Risold, Hallstedt, Gauckler, *Calphad* **20**(3) (1996) 353–361 | SrO(cr), SrO₂(cr), SrO(liq) |
| [3] | Risold, Hallstedt, Gauckler, *J. Am. Ceram. Soc.* **80**(3) (1997) 537–548 | Halite (Ca,Sr)O, liquid L(CaO,SrO) |
| [4] | Jacob, Jayadevan, *J. Phase Equilibria* **21**(4) (2000) 350–356 | Halite (Ca,Sr)O — cross-validation + miscibility gap data |
| [5] | Aljarrah, Medraj, *Calphad* **32** (2008) 240–251 | Metallic Ca-Sr: liquid, fcc, bcc |
| [6] | Dinsdale, *Calphad* **15**(4) (1991) 317–425 | SGTE unaries (copied verbatim from PURE5) |

The original halite data in [4] are from Jacob and Waseda, *J. Am. Ceram. Soc.*
**81**(4) (1998) 1065–1068 (their ref. [8]); [3] fitted those same data.

---

## 3. The problem to be solved

The two starting papers use mutually incompatible liquid models:

- **Ca-O [1]**: associate model (Ca, CaO, O) + 3rd generation — CaO(cr) with three
  Einstein temperatures, CaO(liq) with the **two-state** model.
- **Sr-O [2]**: **two-sublattice ionic** liquid `(Sr+2)₂(Va⁻²,O⁻²)₂` + classical 2nd
  generation.

Two complete binaries were also missing: metallic Ca-Sr and the CaO-SrO pseudobinary.

### Decision: everything to 2nd generation in the liquid phase

The two-state model in Thermo-Calc is a property of the **whole phase**, not of a
species. It cannot coexist with a polynomial SrO(liq) inside the same associate phase:
either the whole phase is two-state, or none of it is. And no 3rd generation
description exists for Sr or for SrO, so forcing 3rd generation would have meant
reassessing SrO from scratch with no data to support it.

The conflict concerns **only the liquid phase**. CaO(cr) and CaO₂ are stoichiometric:
their Einstein functions are written literally in the TDB with `LN` and `EXP` and
clash with nothing, so they are kept as they are. The improvement of [1] in Cp and S°
is not lost.

---

## 4. Three findings from the assembly

### 4.1 The Sr-O liquid needs no reoptimization

An ionic liquid `(A+2)ₚ(Va,O⁻²)_q` with a **single cation** and **ideal mixing** is
mathematically identical to an ideal associate solution {Sr, SrO}:

- `y_Va ≡ x_Sr`, `y_O⁻² ≡ x_SrO`
- The ionic formula unit contains 2 species units, so `2R·Σy·ln y` per formula =
  `R·Σx·ln x` per species — the configurational entropy matches term by term.

The conversion is algebraic. And since [2] used no Sr-SrO interaction parameter (they
assumed it ideal by analogy with Ba-O and Ca-O), `L(LIQUID,SR,SRO) = 0` **exactly**.

### 4.2 Erratum in the published SrO parameters

The `T·ln(T)` coefficient of `G_SrO(cr)` is **47.56**, not 47.36. Some digital
transcriptions of the paper incorrectly gave 47.36.

| Coefficient | S°(298.15) | Δ_fH(298.15) |
|---|---|---|
| 47.36 | 52.27 | −592.20 kJ/mol |
| **47.56** | **53.61** | **−592.14 kJ/mol** |
| Paper, Tables 2 and 3 | 53.58 | −592.15 kJ/mol |

Checked against the original PDF. With 47.36 the whole SrO branch would have been
silently displaced.

### 4.3 Cross-validation of the halite

Two independent routes to the same parameter:

- **[3] Table VI**, published CALPHAD parameter for phase 1x0:
  `L = 23000 − 3·T + 1185·(y_Ca − y_Sr)` J/mol
- **[4] Eqs. 5 and 6**, partial excess Gibbs energies in Hardy subregular form,
  converted to Redlich-Kister in this work:
  `L0 = 25977.5 − 5.58·T`, `L1 = 1062.5 − 0.24·T` J/mol

At 1100 K: L0 = 19700 against 19840 J/mol. **0.7 % difference.** Not full independence
—[3] fitted the data of [4]— but it confirms that the Hardy→RK conversion is correct.
That of [3] was adopted, being the assessed CALPHAD parameter.

A stricter check than a single temperature is the critical point of the gap, which
depends on L0 and L1 together and on their T dependence. Evaluated from the converted
parameters:

| | Tc | x(CaO) at the critical point |
|---|---|---|
| Converted from [4] / Jacob & Waseda | 1172.8 K | 0.530 |
| Measured by Jacob & Waseda | 1173 ± 3 K | 0.53 ± 0.01 |

The conversion reproduces the critical point those authors report, to 0.2 K and 0.001
in mole fraction. Since the critical point is not one of the quantities used to obtain
α and β, this is a genuine check on the algebra rather than a restatement of it.

*Derivation of the conversion:* starting from `G^E = x₁x₂(α·x₁ + β·x₂)` with 1 = CaO,
2 = SrO, the partials give `2α − β = 29165 − 6.3T` and `2β − α = 22790 − 4.86T`, whence
`α = 27040 − 5.82T` and `β = 24915 − 5.34T`, and `L0 = (α+β)/2`, `L1 = (α−β)/2`.
Check: at 1100 K and x = 0.5 it gives G^E = 4960 J/mol, exactly the value the paper
reports.

---

## 5. Contents of the database

### Phases

| Phase | Model | Constituents |
|---|---|---|
| `LIQUID` | associates, 1 site | Ca, CaO, O, Sr, SrO |
| `HALITE` | compound energy, (Ca,Sr)₁(O)₁ | continuous solid solution with a gap |
| `CAO2` | stoichiometric | (Ca)₁(O)₂ |
| `SRO2` | stoichiometric | (Sr)₁(O)₂ |
| `FCC_A1` | (Ca,Sr)₁(Va)₁ | continuous solid solution |
| `BCC_A2` | (Ca,Sr)₁(Va)₃ | continuous solid solution |
| `GAS` | ideal, 1 site | **O₂ only** — incomplete on purpose |

### Own Gibbs energy functions

```
GCAOCR   (298.15–3222)  3T Einstein from [1] Eq. 8
         (3222–6000)    Cp damping from [1] Eq. 10
GCAOLIQ  (298.15–6000)  −656748.5662 + 572.751716·T − 84.370735·T·ln(T)
                        + 2.374069e−4·T²          ← refit in this work
GCAO2CR  (298.15–6000)  3T Einstein from [1] Eq. 8, without the a,b,c terms
GSROCR   (298.15–3000)  −607870 + 268.9·T − 47.56·T·ln(T) − 0.00307·T² + 190000/T
GSROLIQ  (298.15–6000)  −566346 + 449.0·T − 73.1·T·ln(T)
GSRO2CR  (298.15–3000)  GSROCR + GHSEROO − 43740 + 70·T
```

The unaries `GHSERCA`, `GLIQCA`, `GBCCCA`, `GHSERSR`, `GLIQSR`, `GBCCSR`, `GHSEROO`
and `GLIQOO` are copied verbatim from SGTE PURE5.

The Einstein terms are written with the **`GEIN`** function built into Thermo-Calc,
documented in the Database Manager User Guide:

```
GEIN(theta) = 1.5*R*theta + 3*R*T*LN(1-EXP(-theta/T))
```

so that `GCAOCR = E0 + Σ αᵢ·GEIN(θᵢ) − (a/2)·T² − EXP(b+c·T)/c²`. This is the form the
software intends for 3rd generation descriptions, it keeps the lines short, and it
uses the `R = 8.31451` of Thermo-Calc instead of the `8.31446` of [1] — a 6 ppm
difference, negligible (< 0.1 J/mol). The file includes, as a comment, the explicit
form with `LN(1-EXP(...))` in case `GEIN` were unavailable.

### Interaction parameters

| Parameter | Value | Source |
|---|---|---|
| `L(LIQUID,CA,CAO;0)` | −9625 + 12.30·T | [1], accepted version (Fischbach) |
| `L(LIQUID,CAO,O;0)` | +30·T | [1] |
| `L(LIQUID,SR,SRO;0)` | 0 | [2], exact |
| `L(LIQUID,SRO,O;0)` | +30·T | estimated by analogy |
| `L(LIQUID,CAO,SRO;0)` | +25000 | [3], converted |
| `L(LIQUID,CA,SR;0,1,2)` | 1682+0.59T / 521.05+0.052T / −900 | [5] |
| `L(HALITE,CA,SR:O;0)` | 23000 − 3·T | [3] |
| `L(HALITE,CA,SR:O;1)` | +1185 | [3] |
| `L(FCC_A1,CA,SR:VA;0)` | 3770.03 + 0.11·T | [5] |
| `L(BCC_A2,CA,SR:VA;0)` | 3770.03 + 0.01·T | [5] |

From [5] the **random solution model column** (Redlich-Kister) was used, not the
modified quasichemical one, because that one goes straight into an associate TDB.

---

## 6. The CaO liquid refit

`G_CaO(liq)` of [1] is in the two-state model:

```
G_liq-am = G_am − RT·ln(1 + exp(−ΔG_d/RT)),   ΔG_d = A + BT + C·T·ln(T)
```

It was refitted to classical SGTE form `A + B·T + C·T·ln(T) + D·T²` imposing as
**exact constraints**:

- `G_liq(3222) = G_cryst(3222)` → melting point
- `S_liq(Tm) − S_cryst(Tm) = ΔS_fus` → enthalpy of fusion

and fitting the two remaining degrees of freedom by least squares against the
two-state curve over 1000–5000 K.

Result: rms 186 J/mol, maximum deviation 986 J/mol (at the low end, the glass
transition region). `Cp(liq) = 84.371 − 4.748e−4·T`, that is 82.8 J/mol·K at Tm against
the Gurvich estimate of 84 ± 10.

**What is lost:** the glass transition. *Metastable* liquid CaO below ~1500 K no
longer reproduces the rise in Cp on cooling. Irrelevant for slag equilibria; relevant
if one wanted to model vitrification.

**What was verified:** `G_liq − G_cryst > 0` throughout 298–3222 K, that is, the liquid
does not spuriously restabilize at low temperature.

---

## 7. Validation

Done in Python by evaluating the functions **exactly as written in the TDB**.

| Check | Calculated | Reference |
|---|---|---|
| CaO(cr): S°, Cp, Δ_fH | 40.354 / 42.729 / −634.601 | 40.35 / 42.73 / −634.6 |
| CaO₂(cr): S°, Cp, Δ_fH | 59.600 / 61.627 / −648.103 | 59.60 / 61.63 / −648.1 |
| SrO(cr): S°, Δ_fH | 53.606 / −592.143 | 53.58 / −592.15 |
| SrO₂(cr): Δ_fH | −635.883 | −636 |
| Tm(CaO), ΔH_fus | 3222.07 K, 80.923 kJ/mol | 3222 K, 80.92 |
| Tm(SrO), ΔH_fus | 2870.45 K, 89.408 kJ/mol | 2870 K, 89.41 |
| Tf(Ca) pure bcc | 1115.0 K | 1115 K (SGTE) |
| **bcc-Ca + CaO eutectic** | **1093.6 K, 1.78 at.%O** | 1094.6 K, 1.70 at.%O [1] |
| **Halite gap at 1100 K** | **x_SrO = 0.232 / 0.696** | 0.240 / 0.710 experimental [4] |
| Closure of the gap | ~1180 K | ~1170–1220 K |
| G_liq > G_cryst below Tm | no exceptions (CaO and SrO) | — |

The two rows in bold are the ones that matter: they are invariants that were **not**
imposed in the fit, they come out of the complete assembly. The eutectic falls within
1 K and 0.08 at.% of the value of [1]. The miscibility gap reproduces the experimental
data of [4] without having been fitted to them.

### Syntax: rules checked against the documentation

The rules were taken from the **Database Manager's Guide** (Thermo-Calc
Documentation Set 2025a, `Documentation/` or `Thermo-Calc_Documentation-Set.pdf`) and cross-checked
against `PGEO.TDB`, a real database shipped with the software.

| Rule | Source | Status |
|---|---|---|
| Maximum **78 characters** per line | Guide, keyword `FUNCTION` | longest line: 76 |
| Valid unary functions: `LN`, `LOG`, `EXP`, `GEIN` | Guide, `ENTER_PARAMETER` | only those are used |
| `**` only with integer powers | Guide, `ENTER_PARAMETER` | `T**2`, `T**3`, `T**(-1)` |
| Valid phase types: `G A Y L I F B` | Guide, `CONSTITUENT` | `LIQUID:L`, `GAS:G` |
| The type suffix must be repeated in `CONSTITUENT` | Guide, `CONSTITUENT` | correct |
| All `ELEMENT` contiguous, `LIST_OF_REFERENCES` at the end | `PGEO.TDB` | correct |
| ASCII only | — | correct |

A lint script additionally checks that every function invoked with `#` is defined,
that every `REFn` used is declared, that each phase declares as many sites as
sublattices, and that `CONSTITUENT` lists as many groups as sublattices. **It passes
with no errors.**

### Cross-check between the TDB and the macro

| Check | Result |
|---|---|
| `def-sys` of the TCM vs elements of the TDB | CA, O, SR — they match |
| Phases cited in the TCM vs phases of the TDB | all 7, none invented |
| `switch-database` without quotes or path separators | correct |
| `x(...)` conditions on real components | correct |
| Every `MAP` preceded by `SAVE` | 3 of 3 |
| Valid phase statuses (`sus`/`ent`/`fix`/`dor`) | correct |
| Every `=fix 0` releases a condition with `s-c t=none` | 2 of 2 |

### Status in Thermo-Calc

**The database loads and calculates in Thermo-Calc 2025a.** The parser accepts it and
the equilibria in the table above are reproduced from within the software itself, not
only by evaluating the functions in Python.

To reproduce it: set the working directory to this folder and, in Console Mode, run
`macro_file_open validar_CaSrO.TCM`. The 7 blocks are separated by `@&` so they can be
run one at a time.

**It runs in two configurations, from the same file.** There are not two databases: the
difference is which phases are left active with `change-status`.

| Configuration | Phases | What for |
|---|---|---|
| Complete | all 7 | Ca-O diagrams and stability of the peroxides against pO₂. The GAS phase is suspended for the Ca-CaO liquidus at 1 bar (see section 8). |
| CaO-SrO pseudobinary | `LIQUID` and `HALITE` only | `change-status phase BCC_A2 FCC_A1 CAO2 SRO2 GAS = suspended`, section at `X(O) = 0.5`. This is the one in `diagrama_CaSrO_opt.TCM`. |

The peroxides are suspended in the pseudobinary section because at `X(O) = 0.5` they
do not take part and only hinder convergence. The halite gap additionally requires
starting inside it (1100 K, `X(SR) = 0.2`) so that global minimization generates the
two composition sets.

> **Beware of the space in the folder name.** Console Mode `SWITCH_DATABASE` does not
> accept quoted paths, so the macro calls the file by its bare name
> (`sw user CaSrO.TDB`) and the working directory must be set beforehand from
> *File → Set Working Directory*. To avoid this, copy the two files to a folder
> without spaces.

---

## 8. Limitations and warnings

**`L(LIQUID,CAO,SRO) = 25000` is not an experimental datum.** It is the 50000 of [3]
Table VI (given per mole of formula unit of the ionic liquid `(Sr,Ca)₂O₂`, which is two
oxide units) converted to associates. The authors of [3] themselves state that the
mixing behaviour of the SrO-CaO liquid is **unknown** and used that value as a fitting
parameter for melting temperatures in the Sr-Ca-Cu-O quaternary. It is the weakest
piece of the database. The implied critical T (~1500 K) lies far below the CaO-SrO
liquidus (2870–3222 K), so it does not generate a spurious gap in the stable region.

**`L(LIQUID,SRO,O) = +30·T` is an analogy**, copied from the Ca-O value of [1] (which
in turn comes from Liang et al., where it was introduced to fit the decomposition of
CaO on the oxygen-rich side to the measurement of Manara). There is no experimental
datum in Sr-O; in fact [2] **does not extend the liquid beyond the SrO composition**.

**There are no cross parameters** `L(LIQUID,CA,SRO)` or `L(LIQUID,SR,CAO)`. No
information exists and they are left at zero by default.

**The GAS phase is incomplete on purpose.** O₂ only, which suffices to fix pO₂ and to
compute the stability limits of CaO₂ and SrO₂ — in the oxygen-rich region where
decomposition occurs, [1] verifies that O₂ dominates and that any other species stays
below 1e−29 atm. Ca(g), Ca₂(g), O(g), O₃(g), Sr(g), Sr₂(g), Sr₂O(g), CaO(g) and SrO(g)
are missing; they come from SSUB / Lamoreaux. **For the Ca-CaO liquidus at 1 bar the
GAS phase must be suspended**, just as [1] does in its Fig. 12.

**The metals do not dissolve oxygen**, consistent with [1] and [2].

**Effective range of validity**: 298.15 K up to 3000 K (limited by the Sr and SrO
functions of SGTE and [2]); the Ca-O functions reach 6000 K.

---

## 9. Possible improvements

1. **Mixing data for the CaO-SrO liquid.** This is the real gap. Any measurement of
   enthalpy of mixing or of liquidus in the pseudobinary would allow the estimated
   25000 to be replaced by a fitted parameter.
2. **Complete the GAS phase** with SSUB if vaporization equilibria or the Ca-O diagram
   at 1 bar without suspending the gas are needed.
3. **Alternative version of Ca-O**: [1] gives a second set of parameters
   `L(Ca,CaO) = +9284` fitted to the Zaitsev activities instead of the Fischbach
   solubilities. It is noted in the TDB in case a comparison is wanted; it gives a
   eutectic at 1103.2 K and 0.98 at.%O.
4. If the database is to be coupled to an existing slag database, check that the set of
   associates and the reference states match that database exactly.

The primary source of the halite data, Jacob & Waseda (1998), was obtained after the
first version of this database was assembled; its own assessment is compared against
the conversion in section 4.3.

---

## 10. License and how to cite

**Author:** Fernando Ivan Cruz-Cervantes, Instituto Politécnico Nacional.

| What | License |
|---|---|
| Database, data, figures, this README | [CC BY 4.0](LICENSE-DATA.md) |
| Code (`*.py`, `*.TCM`, `*.POP`) | [MIT](LICENSE) |

Citation metadata is in [`CITATION.cff`](CITATION.cff); GitHub generates the *Cite this
repository* button from it.

**If you use this database, cite the original assessments too**, from section 2. What
is original to this work is the assembly: the conversion of the Sr-O ionic liquid to
associates, the refit of the two-state CaO liquid to SGTE form, the Hardy→Redlich-Kister
conversion of the halite, and the parameters marked `REF7` in the TDB. The remaining
parameters belong to their authors.

The unary functions are SGTE data (Dinsdale 1991) and are not covered by the license of
this repository; see the notice in [`LICENSE-DATA.md`](LICENSE-DATA.md).

### Disclaimer

This database is published as is, without warranty. It is academic work, not a
validated commercial database: read section 8 before using it for anything that
matters. In particular, `L(LIQUID,CAO,SRO)` is a fitting parameter, not a measured
datum.
