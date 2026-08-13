# Datos Experimentales Recopilados para Optimización con PARROT (Thermo-Calc)

Este documento centraliza todos los puntos experimentales extraídos de la literatura del sistema **Ca-Sr-O** (`j_calphad_2020`, `risold1996`, `risold1997`, `jacob2000`), estructurados formalmente para facilitar su conversión directa a archivos de entrada `.POP` del módulo **PARROT** de Thermo-Calc.

---

## 1. Fase HALITE $(Ca,Sr)O_1$ - Laguna de Miscibilidad y Líneas de Enlace

### 1.1 Límites de la Laguna de Miscibilidad a 1100 K
* **Referencia:** Jacob (2000) / Jacob & Waseda (1998)
* **Condiciones:** $P = 101325\text{ Pa}$ ($1\text{ bar}$), $T = 1100\text{ K}$, $x(O) = 0.5$
* **Fase:** `HALITE` (Estructura NaCl, $(Ca,Sr)_1O_1$)

| Punto | Componente / Variable | Valor Experimental | Descripción |
|---|---|---|---|
| 1 | $x(SrO)_{\text{límite Ca}}$ | 0.240 | Límite de solubilidad rico en CaO |
| 2 | $x(SrO)_{\text{límite Sr}}$ | 0.712 | Límite de solubilidad rico en SrO |

*Nota en POP:*
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

### 1.2 Líneas de Enlace (Tie-lines) a 1100 K: $(Ca_{1-x}Sr_x)O$ + $(Ca_{1-y}Sr_y)_2PbO_4$
* **Referencia:** Jacob (2000), Tabla 1
* **Condiciones:** $T = 1100\text{ K}$, $P = 101325\text{ Pa}$
* **Variables:** $x = x(SrO)$ en Halita, $y = x(Sr_2PbO_4)$ en la fase ortorrómbica

| ID | $x(SrO)$ en Halita $(Ca,Sr)O$ | $y(Sr_2PbO_4)$ en Ortorrómbica | Estado / Región |
|---|---|---|---|
| TL-01 | 0.006 | 0.088 | Solución sólida continua |
| TL-02 | 0.013 | 0.344 | Solución sólida continua |
| TL-03 | 0.015 | 0.504 | Solución sólida continua |
| TL-04 | 0.018 | 0.653 | Solución sólida continua |
| TL-05 | 0.024 | 0.773 | Solución sólida continua |
| TL-06 | 0.046 | 0.889 | Solución sólida continua |
| TL-07 | 0.073 | 0.926 | Solución sólida continua |
| TL-08 | 0.121 | 0.948 | Solución sólida continua |
| TL-09 | 0.174 | 0.957 | Solución sólida continua |
| TL-10 | 0.240 | 0.961 | **Borde Laguna Miscibilidad (Rico en Ca)** |
| TL-11 | 0.712 | 0.961 | **Borde Laguna Miscibilidad (Rico en Sr)** |
| TL-12 | 0.805 | 0.966 | Solución sólida continua |
| TL-13 | 0.914 | 0.980 | Solución sólida continua |

---

### 1.3 Variación Térmica de la Laguna de Miscibilidad $(Ca,Sr)O$
* **Referencia:** Roth (1989) / Risold et al. (1997), Figura 1
* **Atmósfera:** Aire ($P = 1\text{ bar}$)

| Temperatura $T$ (K) | Fracción Molar SrO ($x_{Ca-rich}$) | Fracción Molar SrO ($x_{Sr-rich}$) | Observación |
|---|---|---|---|
| 873 | ~0.08 | ~0.88 | Medición XRD |
| 973 | ~0.12 | ~0.84 | Medición XRD |
| 1073 | ~0.18 | ~0.78 | Medición XRD |
| 1123 | 0.240 | 0.712 | Medición Jacob (1998, 2000) |
| 1173 | ~0.30 | ~0.65 | Medición XRD |
| 1223 | ~0.40 | ~0.55 | Máximo/Cierre según Roth |
| ~1180–1220 | — | — | **Temperatura consolute (cierre calculado/estimado)** |

---

### 1.4 Entalpía de Mezcla de la Halita $\Delta H_{mix}((Ca_{1-x}Sr_x)O)$
* **Referencia:** Flidlider et al. (1966) / Risold et al. (1997), Figura 3
* **Método:** Calorimetría de disolución en $HClO_4$ a 298.15 K

| $x(SrO)$ | $\Delta H_{mix}$ experimental (kJ/mol) | $\Delta H_{mix}$ calculado TDB (kJ/mol) |
|---|---|---|
| 0.10 | ~1.8 | 1.83 |
| 0.25 | ~3.7 | 3.79 |
| 0.50 | ~5.0 | 4.96 |
| 0.75 | ~3.6 | 3.61 |
| 0.90 | ~1.6 | 1.62 |

---

## 2. Fase LIQUID y Equilibrios en el Binario Ca–O

### 2.1 Solubilidad de CaO en Ca Líquido (Curva Liquidus Ca–CaO)
* **Referencia 1:** Fischbach (1985) (Análisis Químico y ATD)
* **Referencia 2:** Zaitsev et al. (1998, 1999) (Análisis Químico e Inducción de Plasma)

| Temperatura $T$ (K) | Solubilidad $x(CaO)$ (mol%) | Fracción Atómica $x(O)$ | Fuente / Método |
|---|---|---|---|
| 1165 | ~1.0 | ~0.0099 | Zaitsev (1998) |
| 1233 | ~1.2 | ~0.0119 | Fischbach (1985) ATD |
| 1350 | ~2.1 | ~0.0206 | Fischbach (1985) Química |
| 1450 | ~3.2 | ~0.0310 | Fischbach (1985) / Zaitsev |
| 1523 | 16.5 (Descartado: impresciso) | 0.142 | Bevan et al. (1956) (Estimación visual) |
| 1550 | ~4.8 | ~0.0458 | Fischbach (1985) Química |
| 1650 | ~7.2 | ~0.0672 | Fischbach (1985) ATD |
| 1705 | ~9.1 | ~0.0834 | Fischbach (1985) ATD |
| 1723 | ~10.5 | ~0.0950 | Zaitsev (1998) |

---

### 2.2 Actividad de Ca en el Líquido Saturado con CaO
* **Referencia:** Zaitsev et al. (1998, 1999)
* **Método:** Presión de vapor estática / Condensación de Ca

| Temperatura $T$ (K) | Actividad de Ca $a(Ca)_{\text{liq}}$ (Estándar Ca líquido puro) |
|---|---|
| 1165 | 0.992 |
| 1273 | 0.985 |
| 1373 | 0.976 |
| 1473 | 0.963 |
| 1573 | 0.948 |
| 1673 | 0.931 |
| 1723 | 0.920 |

---

### 2.3 Equilibrios Invariantes y Puntos Críticos en Ca–O
* **Eutéctico $\beta\text{-Ca} + CaO \rightleftharpoons L$:**
  * Temperatura: $1107 \pm 1.5\text{ K}$ (Bevan et al. 1956) / $1093.6\text{ K}$ (Calculado TDB).
  * Composición líquida: $x(O) \approx 0.017$ ($1.7\text{ at.\% O}$).
* **Fusión Congruente de CaO(cr):**
  * Temperatura de fusión: $3222 \pm 25\text{ K}$ (Manara et al. 2005 en aire a $0.3\text{ MPa}$).
  * Entalpía de fusión: $\Delta H_{fus} = 80.92\text{ kJ/mol}$ (Alvares et al. 2018 / Deffrennes et al. 2020).

---

## 3. Propiedades Termodinámicas y Puntos Experimentales del Binario Sr–O

### 3.1 Puntos de Fusión y Propiedades Fundamentales de SrO y SrO₂
* **Referencia:** Risold, Hallstedt, Gauckler (1996)

| Propiedad | Valor Experimental / Adoptado | Fuente / Referencia |
|---|---|---|
| $T_m(SrO)$ | $2703\text{ K}$ | Schumacher (1926) |
| $T_m(SrO)$ | $2805\text{ K}$ | Irgashov et al. (1985) |
| $T_m(SrO)$ | $2872\text{ K}$ | Noguchi (1965) |
| $T_m(SrO)$ | $2938\text{ K}$ | Foex (1965) |
| **$T_m(SrO)$ Adoptado** | **$2870\text{ K}$** | **Risold et al. (1996)** |
| $\Delta H_{fus}(SrO)$ | $80.95\text{ kJ/mol}$ (Irgashov) / **$89.41\text{ kJ/mol}$** (Adoptado) | Irgashov (1985) / Risold (1996) |
| $S^\circ_{298}(SrO)$ | $53.63\text{ J/(mol K)}$ | Cordfunke et al. (1994) |
| $\Delta_f H^\circ_{298}(SrO)$ | $-592.15\text{ kJ/mol}$ | Cordfunke et al. (1990) |
| $\Delta_f H^\circ_{298}(SrO_2)$ | $-636\text{ kJ/mol}$ | de Forcrand (1908) / Vedeneev (1952) |

---

### 3.2 Energía Libre de Formación de SrO(s)
* **Referencia:** Ono et al. (1993)
* **Reacción:** $Sr(\text{fcc}) + \frac{1}{2}O_2(\text{g}) \rightleftharpoons SrO(\text{s})$

| Temperatura $T$ (K) | $\Delta_f G^\circ(SrO)$ (kJ/mol) |
|---|---|
| 1373 | -479.2 |
| 1473 | -468.5 |
| 1573 | -457.9 |
| 1673 | -447.2 |
| 1773 | -436.5 |

---

### 3.3 Presión de Descomposición del Peróxido $SrO_2(s) \rightleftharpoons SrO(s) + \frac{1}{2}O_2(g)$
* **Referencia:** Holtermann (1940), Blumenthal (1934, 1935)

| Temperatura $T$ (K) | Presión de Oxígeno $P(O_2)$ (bar) | Estado / Observación |
|---|---|---|
| 488 | ~1.0 | Temperatura de descomposición a 1 bar |
| 673 | 100.0 | Síntesis a alta presión (Rendimiento 16%) |

---

## 4. Resumen de Parámetros Sin Datos Experimentales Directos (Candidatos a Evaluación)

Para conocimiento del módulo/agente de PARROT, los siguientes parámetros **NO disponen de datos experimentales directos** y son estimaciones/analogías en la base actual:

1. `L(LIQUID,CAO,SRO;0) = +25000`:
   * **Origen:** Convertido de Risold (1997) ($+50000\text{ J/mol}$ en 2 subredes iónicas).
   * **Estatus:** Parámetro arbitrario ajustado para reproducir el liquidus en el cuaternario Sr-Ca-Cu-O. El líquido pseudobinario CaO-SrO es experimentalmente **desconocido**.
2. `L(LIQUID,SRO,O;0) = +30*T`:
   * **Origen:** Analogía con Ca-O.
   * **Estatus:** Estimado. No existen datos en el líquido Sr-O ricos en oxígeno.
3. `L(LIQUID,CA,SRO)` y `L(LIQUID,SR,CAO)`:
   * **Estatus:** Faltan parámetros cruzados (asumidos en 0).

---

## 5. Instrucciones de Parseo para el Agente Especializado PARROT

Para generar los archivos `.POP` de PARROT a partir de esta información:
1. Usar la sintaxis estándar de Thermo-Calc POP:
   * `CREATE_EQUILIBRIUM`
   * `SET_CONDITION T=..., P=101325, N=1, X(...)=...`
   * `EXPERIMENT X(...)=...:ERROR` o `EXPERIMENT H=...:ERROR`
   * `SET_START_VALUE`
2. Ponderación sugerida:
   * Datos de tie-lines de la Halita (Sección 1.2): peso / error relativo $\pm 0.005$ en fracción molar.
   * Datos de solubilidad de Ca-O (Sección 2.1): peso / error relativo $\pm 0.002$ en $x(O)$.
   * Puntos de fusión $T_m$: error $\pm 5\text{ K}$.
