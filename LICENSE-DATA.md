# Licencia de la base de datos y los datos

La base de datos termodinámica y los datos de este repositorio —`CaSrO.TDB`,
`CaSrO_opt.tdb`, `CaSrO_opt_liq.tdb`, los archivos `.csv`, `.exp` y `.EXP`, las
figuras y este `README.md`— se distribuyen bajo

**Creative Commons Attribution 4.0 International (CC BY 4.0)**

Copyright (c) 2026 Fernando Cruz

Eres libre de compartir y adaptar este material, incluso comercialmente, siempre que
des el crédito correspondiente (ver `CITATION.cff`), enlaces a la licencia e indiques
si hiciste cambios.

Texto legal completo: <https://creativecommons.org/licenses/by/4.0/legalcode>
Resumen: <https://creativecommons.org/licenses/by/4.0/>

El código (`*.py`, `*.TCM`, `*.POP`) va bajo licencia MIT; ver `LICENSE`.

---

## Procedencia de los parámetros

Esta base **ensambla** descripciones publicadas por otros autores, además de los
parámetros derivados en este trabajo. La licencia CC BY 4.0 cubre el ensamblado, la
conversión entre modelos y los parámetros propios, no la autoría original de los
parámetros tomados de la literatura. Cada uno está atribuido en el
`LIST_OF_REFERENCES` del TDB y en la sección 2 del README. Si usas esta base, cita
también los assessments originales.

## Aviso sobre las funciones unarias SGTE

Las funciones de los elementos puros (`GHSERCA`, `GLIQCA`, `GBCCCA`, `GHSERSR`,
`GLIQSR`, `GBCCSR`, `GHSEROO`, `GLIQOO`) son datos unarios de SGTE, publicados en

> A.T. Dinsdale, "SGTE data for pure elements", *Calphad* **15**(4) (1991) 317–425.

No son obra del autor de este repositorio y **no están cubiertas por la licencia
CC BY 4.0 de arriba**. Su uso y redistribución se rigen por las condiciones de SGTE
(<https://www.sgte.net>). Se incluyen aquí porque una base CALPHAD no es utilizable
sin ellas, con atribución explícita a la fuente.
