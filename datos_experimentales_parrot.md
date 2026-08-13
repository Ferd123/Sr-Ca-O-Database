# Experimental data compiled for optimization with PARROT (Thermo-Calc)

This document collects the experimental points extracted from the literature of the **Ca-Sr-O** system (`j_calphad_2020`, `risold1996`, `risold1997`, `jacob2000`), structured so that they can be converted directly into `.POP` input files for the **PARROT** module of Thermo-Calc.

> **Warning.** Some values in this document are wrong and are **not** used by
> `CaSrO.POP` or by `experimental_CaSrO.EXP`. Specifically, the Ca-CaO liquidus of
> section 2.1 gives 0.99 at.%O at 1165 K, below the eutectic composition
> (1.70 at.%O at 1094.6 K), which is thermodynamically impossible: it appears that the
> discarded Zaitsev points were mixed in under the Fischbach label. The Δ_fG(SrO)
> values of section 3.2 are some 26 kJ/mol too negative. The file is kept as a record
> of the compilation; the data actually used are those of `CaSrO.POP`.

---

## 1. HALITE phase $(Ca,Sr)O_1$ — miscibility gap and tie-lines

### 1.1 Limits of the miscibility gap at 1100 K
* **Reference:** Jacob (2000) / Jacob & Waseda (1998)
* **Conditions:** $P = 101325\text{ Pa}$ ($1\text{ bar}$), $T = 1100\text{ K}$, $x(O) = 0.5$
* **Phase:** `HALITE` (NaCl structure, $(Ca,Sr)_1O_1$)

| Point | Component / variable | Experimental value | Description |
|---|---|---|---|
| 1 | $x(SrO)_{\text{Ca limit}}$ | 0.240 | CaO-rich solubility limit |
| 2 | $x(SrO)_{\text{Sr limit}}$ | 0.712 | SrO-rich solubility limit |

*Note in POP:*
```pop
CREATE_EQUILIBRIUM 1 1
SET_CONDITION T=1100, P=101325, N=1
CHANGE_STATUS PHASE *=SUSPENDED
CHANGE_STATUS PHASE HALITE=ENTERED 1
CHANGE_STATUS PHASE HALITE#2=ENTERED 1
EXPERIMENT X(HALITE#1,SR)=0.240:0.005
EXPERIMENT X(HALITE#2,SR)=0.712:0.005
```

---

### 1.2 Tie-lines at 1100 K: $(Ca_{1-x}Sr_x)O$ + $(Ca_{1-y}Sr_y)_2PbO_4$
* **Reference:** Jacob (2000), Table 1
* **Conditions:** $T = 1100\text{ K}$, $P = 101325\text{ Pa}$
* **Variables:** $x = x(SrO)$ in halite, $y = x(Sr_2PbO_4)$ in the orthorhombic phase

| ID | $x(SrO)$ in halite $(Ca,Sr)O$ | $y(Sr_2PbO_4)$ in orthorhombic | State / region |
|---|---|---|---|
| TL-01 | 0.006 | 0.088 | Continuous solid solution |
| TL-02 | 0.013 | 0.344 | Continuous solid solution |
| TL-03 | 0.015 | 0.504 | Continuous solid solution |
| TL-04 | 0.018 | 0.653 | Continuous solid solution |
| TL-05 | 0.024 | 0.773 | Continuous solid solution |
| TL-06 | 0.046 | 0.889 | Continuous solid solution |
| TL-07 | 0.073 | 0.926 | Continuous solid solution |
| TL-08 | 0.121 | 0.948 | Continuous solid solution |
| TL-09 | 0.174 | 0.957 | Continuous solid solution |
| TL-10 | 0.240 | 0.961 | **Miscibility gap edge (Ca-rich)** |
| TL-11 | 0.712 | 0.961 | **Miscibility gap edge (Sr-rich)** |
| TL-12 | 0.805 | 0.966 | Continuous solid solution |
| TL-13 | 0.914 | 0.980 | Continuous solid solution |

---

### 1.3 Thermal variation of the $(Ca,Sr)O$ miscibility gap
* **Reference:** Roth (1989) / Risold et al. (1997), Figure 1
* **Atmosphere:** air ($P = 1\text{ bar}$)

| Temperature $T$ (K) | SrO mole fraction ($x_{Ca-rich}$) | SrO mole fraction ($x_{Sr-rich}$) | Observation |
|---|---|---|---|
| 873 | ~0.08 | ~0.88 | XRD measurement |
| 973 | ~0.12 | ~0.84 | XRD measurement |
| 1073 | ~0.18 | ~0.78 | XRD measurement |
| 1123 | 0.240 | 0.712 | Jacob measurement (1998, 2000) |
| 1173 | ~0.30 | ~0.65 | XRD measurement |
| 1223 | ~0.40 | ~0.55 | Maximum / closure according to Roth |
| ~1180–1220 | — | — | **Consolute temperature (calculated / estimated closure)** |

---

### 1.4 Halite enthalpy of mixing $\Delta H_{mix}((Ca_{1-x}Sr_x)O)$
* **Reference:** Flidlider et al. (1966) / Risold et al. (1997), Figure 3
* **Method:** solution calorimetry in $HClO_4$ at 298.15 K

| $x(SrO)$ | $\Delta H_{mix}$ experimental (kJ/mol) | $\Delta H_{mix}$ calculated from the TDB (kJ/mol) |
|---|---|---|
| 0.10 | ~1.8 | 1.83 |
| 0.25 | ~3.7 | 3.79 |
| 0.50 | ~5.0 | 4.96 |
| 0.75 | ~3.6 | 3.61 |
| 0.90 | ~1.6 | 1.62 |

---

## 2. LIQUID phase and equilibria in the Ca–O binary

### 2.1 Solubility of CaO in liquid Ca (Ca–CaO liquidus)
* **Reference 1:** Fischbach (1985) (chemical analysis and DTA)
* **Reference 2:** Zaitsev et al. (1998, 1999) (chemical analysis and plasma induction)

**These values are wrong; see the warning at the head of this file.**

| Temperature $T$ (K) | Solubility $x(CaO)$ (mol%) | Atomic fraction $x(O)$ | Source / method |
|---|---|---|---|
| 1165 | ~1.0 | ~0.0099 | Zaitsev (1998) |
| 1233 | ~1.2 | ~0.0119 | Fischbach (1985) DTA |
| 1350 | ~2.1 | ~0.0206 | Fischbach (1985) chemical |
| 1450 | ~3.2 | ~0.0310 | Fischbach (1985) / Zaitsev |
| 1523 | 16.5 (discarded: imprecise) | 0.142 | Bevan et al. (1956) (visual estimate) |
| 1550 | ~4.8 | ~0.0458 | Fischbach (1985) chemical |
| 1650 | ~7.2 | ~0.0672 | Fischbach (1985) DTA |
| 1705 | ~9.1 | ~0.0834 | Fischbach (1985) DTA |
| 1723 | ~10.5 | ~0.0950 | Zaitsev (1998) |

---

### 2.2 Activity of Ca in the CaO-saturated liquid
* **Reference:** Zaitsev et al. (1998, 1999)
* **Method:** static vapour pressure / Ca condensation

| Temperature $T$ (K) | Ca activity $a(Ca)_{\text{liq}}$ (standard state: pure liquid Ca) |
|---|---|
| 1165 | 0.992 |
| 1273 | 0.985 |
| 1373 | 0.976 |
| 1473 | 0.963 |
| 1573 | 0.948 |
| 1673 | 0.931 |
| 1723 | 0.920 |

---

### 2.3 Invariant equilibria and critical points in Ca–O
* **Eutectic $\beta\text{-Ca} + CaO \rightleftharpoons L$:**
  * Temperature: $1107 \pm 1.5\text{ K}$ (Bevan et al. 1956) / $1093.6\text{ K}$ (calculated from the TDB).
  * Liquid composition: $x(O) \approx 0.017$ ($1.7\text{ at.\% O}$).
* **Congruent melting of CaO(cr):**
  * Melting temperature: $3222 \pm 25\text{ K}$ (Manara et al. 2005, in air at $0.3\text{ MPa}$).
  * Enthalpy of fusion: $\Delta H_{fus} = 80.92\text{ kJ/mol}$ (Alvares et al. 2018 / Deffrennes et al. 2020).

---

## 3. Thermodynamic properties and experimental points of the Sr–O binary

### 3.1 Melting points and fundamental properties of SrO and SrO₂
* **Reference:** Risold, Hallstedt, Gauckler (1996)

| Property | Experimental / adopted value | Source / reference |
|---|---|---|
| $T_m(SrO)$ | $2703\text{ K}$ | Schumacher (1926) |
| $T_m(SrO)$ | $2805\text{ K}$ | Irgashov et al. (1985) |
| $T_m(SrO)$ | $2872\text{ K}$ | Noguchi (1965) |
| $T_m(SrO)$ | $2938\text{ K}$ | Foex (1965) |
| **$T_m(SrO)$ adopted** | **$2870\text{ K}$** | **Risold et al. (1996)** |
| $\Delta H_{fus}(SrO)$ | $80.95\text{ kJ/mol}$ (Irgashov) / **$89.41\text{ kJ/mol}$** (adopted) | Irgashov (1985) / Risold (1996) |
| $S^\circ_{298}(SrO)$ | $53.63\text{ J/(mol K)}$ | Cordfunke et al. (1994) |
| $\Delta_f H^\circ_{298}(SrO)$ | $-592.15\text{ kJ/mol}$ | Cordfunke et al. (1990) |
| $\Delta_f H^\circ_{298}(SrO_2)$ | $-636\text{ kJ/mol}$ | de Forcrand (1908) / Vedeneev (1952) |

---

### 3.2 Gibbs energy of formation of SrO(s)
* **Reference:** Ono et al. (1993)
* **Reaction:** $Sr(\text{fcc}) + \frac{1}{2}O_2(\text{g}) \rightleftharpoons SrO(\text{s})$

**These values are some 26 kJ/mol too negative; see the warning at the head of this
file.** The database gives −452.6 kJ/mol at 1373 K, against the Risold curve at ~−453.

| Temperature $T$ (K) | $\Delta_f G^\circ(SrO)$ (kJ/mol) |
|---|---|
| 1373 | -479.2 |
| 1473 | -468.5 |
| 1573 | -457.9 |
| 1673 | -447.2 |
| 1773 | -436.5 |

---

### 3.3 Decomposition pressure of the peroxide $SrO_2(s) \rightleftharpoons SrO(s) + \frac{1}{2}O_2(g)$
* **Reference:** Holtermann (1940), Blumenthal (1934, 1935)

| Temperature $T$ (K) | Oxygen pressure $P(O_2)$ (bar) | State / observation |
|---|---|---|
| 488 | ~1.0 | Decomposition temperature at 1 bar |
| 673 | 100.0 | High-pressure synthesis (16% yield) |

---

## 4. Summary of parameters with no direct experimental data (candidates for assessment)

The following parameters **have no direct experimental data** and are estimates or
analogies in the current database:

1. `L(LIQUID,CAO,SRO;0) = +25000`:
   * **Origin:** converted from Risold (1997) ($+50000\text{ J/mol}$ on 2 ionic sublattices).
   * **Status:** arbitrary parameter fitted to reproduce the liquidus in the Sr-Ca-Cu-O quaternary. The CaO-SrO pseudobinary liquid is experimentally **unknown**.
2. `L(LIQUID,SRO,O;0) = +30*T`:
   * **Origin:** analogy with Ca-O.
   * **Status:** estimated. No data exist for the oxygen-rich Sr-O liquid.
3. `L(LIQUID,CA,SRO)` and `L(LIQUID,SR,CAO)`:
   * **Status:** cross parameters missing (assumed zero).

---

## 5. Notes for turning this into a POP file

To generate the PARROT `.POP` files from this information:

1. Use standard Thermo-Calc POP syntax:
   * `CREATE_EQUILIBRIUM`
   * `SET_CONDITION T=..., P=101325, N=1, X(...)=...`
   * `EXPERIMENT X(...)=...:ERROR` or `EXPERIMENT H=...:ERROR`
   * `SET_START_VALUE`
2. Suggested weighting:
   * Halite tie-line data (section 1.2): weight / relative error $\pm 0.005$ in mole fraction.
   * Ca-O solubility data (section 2.1): weight / relative error $\pm 0.002$ in $x(O)$.
   * Melting points $T_m$: error $\pm 5\text{ K}$.

See `CaSrO.POP` for the file actually used, whose construction rules and provenance
notes supersede the suggestions above.
