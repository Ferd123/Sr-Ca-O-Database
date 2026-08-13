# Base de datos termodinámica Ca-Sr-O para Thermo-Calc

Descripción CALPHAD del sistema Ca-Sr-O con la fase líquida en **modelo de asociados**
(especies Ca, Sr, O, CaO, SrO), ensamblada a partir de cinco fuentes publicadas que
usaban modelos mutuamente incompatibles.

**Autor:** Fernando Cruz (Instituto Politécnico Nacional)
Fecha de ensamblado: 2026-08-03. Thermo-Calc 2025a.
Licencia: CC BY 4.0 (base de datos y datos) / MIT (código). Ver sección 10.

---

## 1. Archivos

| Archivo | Qué es |
|---|---|
| `CaSrO.TDB` | La base de datos. Fases condensadas + GAS parcial. |
| `CaSrO_opt.tdb`, `CaSrO_opt_liq.tdb` | Variantes con el líquido reoptimizado. |
| `validar_CaSrO.TCM` | Macro de Console Mode, 7 bloques (ver abajo). |
| `CaSrO.POP`, `optimizar_CaSrO_PARROT.TCM` | Optimización con PARROT. |
| `*.csv`, `*.exp`, `*.EXP` | Datos experimentales de la literatura, con su cita. |
| `*.py` | Scripts de evaluación y de las figuras. |

Los artículos de los que salen los parámetros **no** están en el repositorio: tienen
copyright. Las citas completas están en la sección 2 y en el `LIST_OF_REFERENCES` del
propio TDB.

Los bloques del macro, separados por `@&` para poder correrlos de uno en uno:

| # | Qué hace | Qué debe salir |
|---|---|---|
| 1 | Carga la base y lista el sistema | 7 fases |
| 2 | Fusión congruente de CaO | 3222.0 K |
| 3 | Fusión congruente de SrO | 2870.5 K |
| 4 | Diagrama Ca-O (GAS suspendida) | eutéctico ~1094 K, 1.7 at.%O |
| 5 | Pseudobinario CaO-SrO | laguna x(Sr) 0.116–0.348 a 1100 K |
| 6 | Binario metálico Ca-Sr | fcc y bcc continuas, sin fases intermedias |
| 7 | Estabilidad de CaO₂ frente a pO₂ | descomposición a baja T |

---

## 2. Fuentes

| Ref | Cita | Qué aporta |
|---|---|---|
| [1] | Deffrennes, Jakse, Alvares, Nuta, Pasturel, Khvan, Pisch, *Calphad* **69** (2020) 101764 | CaO(cr), CaO₂(cr), CaO(liq), L(Ca,CaO), L(CaO,O) |
| [2] | Risold, Hallstedt, Gauckler, *Calphad* **20**(3) (1996) 353–361 | SrO(cr), SrO₂(cr), SrO(liq) |
| [3] | Risold, Hallstedt, Gauckler, *J. Am. Ceram. Soc.* **80**(3) (1997) 537–548 | Halita (Ca,Sr)O, L(CaO,SrO) líquido |
| [4] | Jacob, Jayadevan, *J. Phase Equilibria* **21**(4) (2000) 350–356 | Halita (Ca,Sr)O — validación cruzada + datos de la laguna |
| [5] | Aljarrah, Medraj, *Calphad* **32** (2008) 240–251 | Ca-Sr metálico: líquido, fcc, bcc |
| [6] | Dinsdale, *Calphad* **15**(4) (1991) 317–425 | Unarias SGTE (copiadas literalmente de PURE5) |

Los datos originales de la halita en [4] son de Jacob y Waseda, *J. Am. Ceram. Soc.*
**81** (1998) 1065 (su ref. [8]); [3] ajustó a esos mismos datos.

---

## 3. El problema que había que resolver

Los dos papers de partida usan modelos de líquido incompatibles entre sí:

- **Ca-O [1]**: modelo de asociados (Ca, CaO, O) + 3ª generación — CaO(cr) con tres
  temperaturas de Einstein, CaO(liq) con modelo de **dos estados** (two-state).
- **Sr-O [2]**: líquido **iónico de dos subredes** `(Sr+2)₂(Va⁻²,O⁻²)₂` + 2ª generación
  clásica.

Además faltaban dos binarios completos: Ca-Sr metálico y el pseudobinario CaO-SrO.

### Decisión: todo a 2ª generación en la fase líquida

El modelo two-state en Thermo-Calc es una propiedad de la **fase completa**, no de una
especie. No puede coexistir con SrO(liq) polinomial dentro de la misma fase de
asociados: o toda la fase es two-state, o ninguna. Y no existe descripción de 3ª
generación para Sr ni para SrO, así que forzar la 3ª generación implicaba reasesar SrO
desde cero sin datos que lo soporten.

El conflicto es **sólo de la fase líquida**. CaO(cr) y CaO₂ son estequiométricas: sus
funciones de Einstein se escriben literales en el TDB con `LN` y `EXP` y no chocan con
nada, así que se conservan tal cual. No se pierde la mejora de [1] en Cp y S°.

---

## 4. Los tres hallazgos del ensamblado

### 4.1 El líquido Sr-O no necesita reoptimización

Un líquido iónico `(A+2)ₚ(Va,O⁻²)_q` con **un solo catión** y **mezcla ideal** es
matemáticamente idéntico a una solución ideal de asociados {Sr, SrO}:

- `y_Va ≡ x_Sr`, `y_O⁻² ≡ x_SrO`
- La unidad de fórmula iónica contiene 2 unidades de especie, así que
  `2R·Σy·ln y` por fórmula = `R·Σx·ln x` por especie — la entropía configuracional
  coincide término a término.

La conversión es algebraica. Y como [2] no usó ningún parámetro de interacción Sr-SrO
(lo asumieron ideal por analogía con Ba-O y Ca-O), `L(LIQUID,SR,SRO) = 0` **exactamente**.

### 4.2 Errata en los parámetros publicados de SrO

El coeficiente de `T·ln(T)` de `G_SrO(cr)` es **47.56**, no 47.36. El markdown parseado
del paper (OCR de un escaneo de 1996) daba 47.36.

| Coeficiente | S°(298.15) | Δ_fH(298.15) |
|---|---|---|
| 47.36 | 52.27 | −592.20 kJ/mol |
| **47.56** | **53.61** | **−592.14 kJ/mol** |
| Paper, Tablas 2 y 3 | 53.58 | −592.15 kJ/mol |

Verificado contra el PDF original. Con 47.36 toda la rama de SrO habría quedado
desplazada en silencio.

### 4.3 Validación cruzada de la halita

Dos rutas independientes al mismo parámetro:

- **[3] Tabla VI**, parámetro CALPHAD publicado para la fase 1x0:
  `L = 23000 − 3·T + 1185·(y_Ca − y_Sr)` J/mol
- **[4] Ecs. 5 y 6**, energías de Gibbs de exceso parciales en forma subregular de
  Hardy, convertidas a Redlich-Kister en este trabajo:
  `L0 = 25977.5 − 5.58·T`, `L1 = 1062.5 − 0.24·T` J/mol

A 1100 K: L0 = 19700 contra 19840 J/mol. **0.7 % de diferencia.** No es independencia
total —[3] ajustó a los datos de [4]— pero confirma que la conversión Hardy→RK está
bien hecha. Se adoptó la de [3] por ser el parámetro CALPHAD evaluado.

*Derivación de la conversión:* partiendo de `G^E = x₁x₂(α·x₁ + β·x₂)` con 1 = CaO,
2 = SrO, las parciales dan `2α − β = 29165 − 6.3T` y `2β − α = 22790 − 4.86T`, de donde
`α = 27040 − 5.82T` y `β = 24915 − 5.34T`, y `L0 = (α+β)/2`, `L1 = (α−β)/2`.
Comprobación: a 1100 K y x = 0.5 da G^E = 4960 J/mol, exactamente el valor que reporta
el paper.

---

## 5. Contenido de la base

### Fases

| Fase | Modelo | Constituyentes |
|---|---|---|
| `LIQUID` | asociados, 1 sitio | Ca, CaO, O, Sr, SrO |
| `HALITE` | energía de compuesto, (Ca,Sr)₁(O)₁ | solución sólida continua con laguna |
| `CAO2` | estequiométrica | (Ca)₁(O)₂ |
| `SRO2` | estequiométrica | (Sr)₁(O)₂ |
| `FCC_A1` | (Ca,Sr)₁(Va)₁ | solución sólida completa |
| `BCC_A2` | (Ca,Sr)₁(Va)₃ | solución sólida completa |
| `GAS` | ideal, 1 sitio | **sólo O₂** — incompleta a propósito |

### Funciones de Gibbs propias

```
GCAOCR   (298.15–3222)  Einstein 3T de [1] Ec. 8
         (3222–6000)    amortiguamiento de Cp de [1] Ec. 10
GCAOLIQ  (298.15–6000)  −656748.5662 + 572.751716·T − 84.370735·T·ln(T)
                        + 2.374069e−4·T²          ← refit de este trabajo
GCAO2CR  (298.15–6000)  Einstein 3T de [1] Ec. 8, sin términos a,b,c
GSROCR   (298.15–3000)  −607870 + 268.9·T − 47.56·T·ln(T) − 0.00307·T² + 190000/T
GSROLIQ  (298.15–6000)  −566346 + 449.0·T − 73.1·T·ln(T)
GSRO2CR  (298.15–3000)  GSROCR + GHSEROO − 43740 + 70·T
```

Las unarias `GHSERCA`, `GLIQCA`, `GBCCCA`, `GHSERSR`, `GLIQSR`, `GBCCSR`, `GHSEROO`,
`GLIQOO` están copiadas literalmente de SGTE PURE5.

Los términos de Einstein se escriben con la función **`GEIN`** que Thermo-Calc trae
integrada, documentada en el Database Manager User Guide:

```
GEIN(theta) = 1.5*R*theta + 3*R*T*LN(1-EXP(-theta/T))
```

de modo que `GCAOCR = E0 + Σ αᵢ·GEIN(θᵢ) − (a/2)·T² − EXP(b+c·T)/c²`. Es la forma
prevista por el software para descripciones de 3ª generación, mantiene las líneas
cortas y usa el `R = 8.31451` de Thermo-Calc en vez del `8.31446` de [1] — una
diferencia de 6 ppm, inapreciable (< 0.1 J/mol). El archivo incluye, en comentario,
la forma explícita con `LN(1-EXP(...))` por si `GEIN` no estuviera disponible.

### Parámetros de interacción

| Parámetro | Valor | Fuente |
|---|---|---|
| `L(LIQUID,CA,CAO;0)` | −9625 + 12.30·T | [1], versión aceptada (Fischbach) |
| `L(LIQUID,CAO,O;0)` | +30·T | [1] |
| `L(LIQUID,SR,SRO;0)` | 0 | [2], exacto |
| `L(LIQUID,SRO,O;0)` | +30·T | estimado por analogía |
| `L(LIQUID,CAO,SRO;0)` | +25000 | [3], convertido |
| `L(LIQUID,CA,SR;0,1,2)` | 1682+0.59T / 521.05+0.052T / −900 | [5] |
| `L(HALITE,CA,SR:O;0)` | 23000 − 3·T | [3] |
| `L(HALITE,CA,SR:O;1)` | +1185 | [3] |
| `L(FCC_A1,CA,SR:VA;0)` | 3770.03 + 0.11·T | [5] |
| `L(BCC_A2,CA,SR:VA;0)` | 3770.03 + 0.01·T | [5] |

De [5] se usó la **columna del modelo de solución aleatoria** (Redlich-Kister), no la
del cuasiquímico modificado, porque esa sí entra directo en un TDB de asociados.

---

## 6. El refit del líquido CaO

`G_CaO(liq)` de [1] está en modelo de dos estados:

```
G_liq-am = G_am − RT·ln(1 + exp(−ΔG_d/RT)),   ΔG_d = A + BT + C·T·ln(T)
```

Se reajustó a forma SGTE clásica `A + B·T + C·T·ln(T) + D·T²` imponiendo como
**restricciones exactas**:

- `G_liq(3222) = G_cryst(3222)` → punto de fusión
- `S_liq(Tm) − S_cryst(Tm) = ΔS_fus` → entalpía de fusión

y ajustando los dos grados de libertad restantes por mínimos cuadrados contra la curva
two-state en 1000–5000 K.

Resultado: rms 186 J/mol, desviación máxima 986 J/mol (en el extremo bajo, zona de la
transición vítrea). `Cp(liq) = 84.371 − 4.748e−4·T`, o sea 82.8 J/mol·K en Tm contra la
estimación de Gurvich de 84 ± 10.

**Qué se pierde:** la transición vítrea. El CaO líquido *metaestable* por debajo de
~1500 K ya no reproduce el aumento de Cp al enfriar. Irrelevante para equilibrios de
escoria; relevante si alguien quisiera modelar vitrificación.

**Qué se verificó:** `G_liq − G_cryst > 0` en todo 298–3222 K, o sea el líquido no se
reestabiliza espuriamente a baja temperatura.

---

## 7. Validación

Hecha en Python evaluando las funciones **tal como quedaron escritas en el TDB**.

| Comprobación | Calculado | Referencia |
|---|---|---|
| CaO(cr): S°, Cp, Δ_fH | 40.354 / 42.729 / −634.601 | 40.35 / 42.73 / −634.6 |
| CaO₂(cr): S°, Cp, Δ_fH | 59.600 / 61.627 / −648.103 | 59.60 / 61.63 / −648.1 |
| SrO(cr): S°, Δ_fH | 53.606 / −592.143 | 53.58 / −592.15 |
| SrO₂(cr): Δ_fH | −635.883 | −636 |
| Tm(CaO), ΔH_fus | 3222.07 K, 80.923 kJ/mol | 3222 K, 80.92 |
| Tm(SrO), ΔH_fus | 2870.45 K, 89.408 kJ/mol | 2870 K, 89.41 |
| Tf(Ca) bcc puro | 1115.0 K | 1115 K (SGTE) |
| **Eutéctico bcc-Ca + CaO** | **1093.6 K, 1.78 at.%O** | 1094.6 K, 1.70 at.%O [1] |
| **Laguna halita a 1100 K** | **x_SrO = 0.232 / 0.696** | 0.240 / 0.710 experimental [4] |
| Cierre de la laguna | ~1180 K | ~1170–1220 K |
| G_liq > G_cryst bajo Tm | sin excepciones (CaO y SrO) | — |

Las dos filas en negrita son las que importan: son invariantes que **no** se impusieron
en el ajuste, salen del ensamblado completo. El eutéctico cae a 1 K y 0.08 at.% del
valor de [1]. La laguna de miscibilidad reproduce los datos experimentales de [4] sin
haber sido ajustada a ellos.

### Sintaxis: reglas verificadas contra la documentación

Las reglas se sacaron del **Database Manager User Guide** (Thermo-Calc Documentation
Set 2025a, `Manuals/All Thermo-Calc Documentation/`) y se contrastaron con `PGEO.TDB`,
una base real que trae el software.

| Regla | Fuente | Estado |
|---|---|---|
| Máximo **78 caracteres** por línea | Guide, keyword `FUNCTION` | línea más larga: 76 |
| Funciones unarias válidas: `LN`, `LOG`, `EXP`, `GEIN` | Guide, `ENTER_PARAMETER` | sólo se usan ésas |
| `**` únicamente con potencias enteras | Guide, `ENTER_PARAMETER` | `T**2`, `T**3`, `T**(-1)` |
| Tipos de fase válidos: `G A Y L I F B` | Guide, `CONSTITUENT` | `LIQUID:L`, `GAS:G` |
| El sufijo de tipo debe repetirse en `CONSTITUENT` | Guide, `CONSTITUENT` | correcto |
| Todos los `ELEMENT` contiguos, `LIST_OF_REFERENCES` al final | `PGEO.TDB` | correcto |
| Sólo ASCII | — | correcto |

Un script de lint comprueba además que toda función invocada con `#` esté definida,
que toda referencia `REFn` usada esté declarada, que cada fase declare tantos sitios
como subredes, y que `CONSTITUENT` liste tantos grupos como subredes. **Pasa sin
errores.**

### Cruce entre el TDB y el macro

| Comprobación | Resultado |
|---|---|
| `def-sys` del TCM vs elementos del TDB | CA, O, SR — coinciden |
| Fases citadas en el TCM vs fases del TDB | las 7, ninguna inventada |
| `switch-database` sin comillas ni separadores de ruta | correcto |
| Condiciones `x(...)` sobre componentes reales | correcto |
| Cada `MAP` precedido de `SAVE` | 3 de 3 |
| Estados de fase válidos (`sus`/`ent`/`fix`/`dor`) | correctos |
| Cada `=fix 0` libera una condición con `s-c t=none` | 2 de 2 |

### Estado en Thermo-Calc

**La base carga y calcula en Thermo-Calc 2025a.** El parser la acepta y los equilibrios
de la tabla de arriba se reproducen desde el propio software, no sólo evaluando las
funciones en Python.

Para reproducirlo: fijar el directorio de trabajo en esta carpeta y, en Console Mode,
ejecutar `macro_file_open validar_CaSrO.TCM`. Los 7 bloques van separados por `@&`
para poder correrlos de uno en uno.

**Se corre en dos configuraciones, con el mismo archivo.** No hay dos bases: la
diferencia está en qué fases se dejan activas con `change-status`.

| Configuración | Fases | Para qué |
|---|---|---|
| Completa | las 7 | Diagramas Ca-O y estabilidad de los peróxidos frente a pO₂. La fase GAS se suspende para el liquidus Ca-CaO a 1 bar (ver sección 8). |
| Pseudobinario CaO-SrO | sólo `LIQUID` y `HALITE` | `change-status phase BCC_A2 FCC_A1 CAO2 SRO2 GAS = suspended`, corte a `X(O) = 0.5`. Es la de `diagrama_CaSrO_opt.TCM`. |

Los peróxidos se suspenden en el corte pseudobinario porque a `X(O) = 0.5` no
participan y sólo estorban a la convergencia. La laguna de la halita pide además
arrancar dentro de ella (1100 K, `X(SR) = 0.2`) para que la minimización global genere
los dos conjuntos de composición.

> **Ojo con el espacio en el nombre de la carpeta.** `SWITCH_DATABASE` de Console Mode
> no acepta rutas entre comillas, así que el macro llama al archivo por su nombre a
> secas (`sw user CaSrO.TDB`) y hay que fijar antes el directorio de trabajo desde
> *File → Set Working Directory*. Si prefieres evitarlo, copia los dos archivos a una
> carpeta sin espacios.

---

## 8. Limitaciones y advertencias

**`L(LIQUID,CAO,SRO) = 25000` no es un dato experimental.** Es el 50000 de [3]
Tabla VI (dado por mol de unidad de fórmula del líquido iónico `(Sr,Ca)₂O₂`, que son
dos unidades de óxido) convertido a asociados. Los propios autores de [3] declaran que
el comportamiento de mezcla del líquido SrO-CaO es **desconocido** y usaron ese valor
como parámetro de ajuste de temperaturas de fusión en el cuaternario Sr-Ca-Cu-O. Es la
pieza más débil de la base. La T crítica implícita (~1500 K) queda muy por debajo del
liquidus CaO-SrO (2870–3222 K), así que no genera laguna espuria en la región estable.

**`L(LIQUID,SRO,O) = +30·T` es una analogía**, copiada del valor de Ca-O de [1] (que a
su vez viene de Liang et al., donde se introdujo para ajustar la descomposición de CaO
del lado rico en oxígeno a la medición de Manara). No hay dato experimental en Sr-O;
de hecho [2] **no extiende el líquido más allá de la composición SrO**.

**No hay parámetros cruzados** `L(LIQUID,CA,SRO)` ni `L(LIQUID,SR,CAO)`. No existe
información y quedan en cero por defecto.

**La fase GAS está incompleta a propósito.** Sólo O₂, suficiente para fijar pO₂ y
calcular los límites de estabilidad de CaO₂ y SrO₂ — en la región rica en oxígeno donde
ocurre la descomposición, [1] verifica que O₂ domina y que cualquier otra especie está
por debajo de 1e−29 atm. Faltan Ca(g), Ca₂(g), O(g), O₃(g), Sr(g), Sr₂(g), Sr₂O(g),
CaO(g), SrO(g), que vienen de SSUB / Lamoreaux. **Para el liquidus Ca-CaO a 1 bar hay
que suspender la fase GAS**, igual que hace [1] en su Fig. 12.

**Los metales no disuelven oxígeno**, consistente con [1] y [2].

**Rango de validez efectivo**: 298.15 K hasta 3000 K (limitado por las funciones de Sr
y SrO de SGTE y [2]); las funciones de Ca-O llegan a 6000 K.

---

## 9. Mejoras posibles

1. **Conseguir Jacob & Waseda 1998** (*J. Am. Ceram. Soc.* **81**, 1065) para tener la
   evaluación directa de la halita en vez de la cita en [4].
2. **Datos de mezcla del líquido CaO-SrO.** Es el hueco real. Cualquier medición de
   entalpía de mezcla o de liquidus en el pseudobinario permitiría reemplazar el 25000
   estimado por un parámetro ajustado.
3. **Completar la fase GAS** con SSUB si se necesitan equilibrios de vaporización o
   el diagrama Ca-O a 1 bar sin suspender el gas.
4. **Versión alternativa del Ca-O**: [1] da un segundo juego de parámetros
   `L(Ca,CaO) = +9284` ajustado a las actividades de Zaitsev en vez de las
   solubilidades de Fischbach. Está anotado en el TDB por si se quiere comparar; da un
   eutéctico a 1103.2 K y 0.98 at.%O.
5. Si la base se va a acoplar a una de escorias existente, verificar que el conjunto de
   asociados y los estados de referencia calcen exactamente con esa base.

---

## 10. Licencia y cómo citar

**Autor:** Fernando Cruz, Instituto Politécnico Nacional.

| Qué | Licencia |
|---|---|
| Base de datos, datos, figuras, este README | [CC BY 4.0](LICENSE-DATA.md) |
| Código (`*.py`, `*.TCM`, `*.POP`) | [MIT](LICENSE) |

Los metadatos de cita están en [`CITATION.cff`](CITATION.cff); GitHub genera desde ahí
el botón *Cite this repository*.

**Si usas esta base, cita también los assessments originales** de la sección 2. Lo que
es original de este trabajo es el ensamblado: la conversión del líquido iónico Sr-O a
asociados, el refit del líquido CaO de dos estados a forma SGTE, la conversión
Hardy→Redlich-Kister de la halita, y los parámetros marcados `REF7` en el TDB. Los
demás parámetros son de sus autores.

Las funciones unarias son datos de SGTE (Dinsdale 1991) y no están cubiertas por la
licencia de este repositorio; ver el aviso en [`LICENSE-DATA.md`](LICENSE-DATA.md).

### Descargo

Esta base se publica tal cual, sin garantía. Es trabajo académico, no una base
comercial validada: lee la sección 8 antes de usarla para nada que importe. En
particular, `L(LIQUID,CAO,SRO)` es un parámetro de ajuste, no un dato medido.
