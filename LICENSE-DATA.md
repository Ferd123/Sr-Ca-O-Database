# License of the database and the data

The thermodynamic database and the data in this repository —`CaSrO.TDB`,
`CaSrO_opt.tdb`, `CaSrO_opt_liq.tdb`, the `.csv`, `.exp` and `.EXP` files, the figures
and this `README.md`— are distributed under

**Creative Commons Attribution 4.0 International (CC BY 4.0)**

Copyright (c) 2026 Fernando Ivan Cruz-Cervantes

You are free to share and adapt this material, including commercially, as long as you
give appropriate credit (see `CITATION.cff`), link to the license and indicate whether
you made changes.

Full legal text: <https://creativecommons.org/licenses/by/4.0/legalcode>
Summary: <https://creativecommons.org/licenses/by/4.0/>

The macros (`*.TCM`, `*.POP`) are under the MIT license; see `LICENSE`.

---

## Provenance of the parameters

This database **assembles** descriptions published by other authors, in addition to the
parameters derived in this work. The CC BY 4.0 license covers the assembly, the
conversion between models and the parameters of this work, not the original authorship
of the parameters taken from the literature. Each one is attributed in the
`LIST_OF_REFERENCES` of the TDB and in section 2 of the README. If you use this
database, cite the original assessments too.

## Notice on the SGTE unary functions

The pure element functions (`GHSERCA`, `GLIQCA`, `GBCCCA`, `GHSERSR`, `GLIQSR`,
`GBCCSR`, `GHSEROO`, `GLIQOO`) are SGTE unary data, published in

> A.T. Dinsdale, "SGTE data for pure elements", *Calphad* **15**(4) (1991) 317–425.

They are not the work of the author of this repository and are **not covered by the
CC BY 4.0 license above**. Their use and redistribution are governed by the SGTE
conditions (<https://www.sgte.net>). They are included here because a CALPHAD database
is unusable without them, with explicit attribution to the source.
