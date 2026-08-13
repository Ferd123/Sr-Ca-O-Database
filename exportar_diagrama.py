r"""
Converts a DATAPLOT file (.exp) exported by Thermo-Calc into a tidy CSV,
suitable for plotting with any tool.

    python exportar_diagrama.py dos.exp diagrama_CaSrO.csv

The .exp holds the already calculated diagram as numerical coordinates,
grouped into blocks. Each block is one continuous line of the diagram and is
preceded by comments declaring the phases involved:

    $ BLOCK #68 1 FOR:
    $E LIQUID          <- phase that appears or disappears
    $F0 HALITE#2       <- phase fixed at zero amount (the boundary)

The resulting CSV carries one row per point and a column identifying which
boundary it belongs to, so that whoever receives it can draw each line
separately without having Thermo-Calc.
"""
import csv
import io
import re
import sys

ENT = sys.argv[1] if len(sys.argv) > 1 else "dos.exp"
SAL = sys.argv[2] if len(sys.argv) > 2 else "diagrama_CaSrO.csv"

lineas = io.open(ENT, encoding="latin-1").read().splitlines()

meta = {}
for l in lineas:                       # global header
    m = re.match(r"\s*(XTEXT|YTEXT|TITLE)\s+(.*)", l)
    if m:
        meta[m.group(1)] = m.group(2).strip().rstrip(",")

filas = []
bloque = 0
fases = []
dentro = False

for l in lineas:
    s = l.strip()
    if s.startswith("$ BLOCK"):
        bloque += 1
        fases = []
        continue
    m = re.match(r"\$([EF]\d*)\s+(\S+)", s)          # $E / $F0 declare phases
    if m:
        fases.append(f"{m.group(1)}:{m.group(2)}")
        continue
    if s.startswith("BLOCK "):
        dentro = True
        continue
    if s.startswith("BLOCKEND"):
        dentro = False
        continue
    if dentro:
        p = s.split()
        if len(p) >= 2:
            try:
                x, y = float(p[0]), float(p[1])
            except ValueError:
                continue
            filas.append([bloque, " + ".join(fases) or "?", x, y])

with io.open(SAL, "w", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    # Comment lines must not contain commas: the csv module would quote them
    # and they would stop being recognised as comments when the file is read
    # with pandas or a spreadsheet.
    w.writerow([f"# {meta.get('TITLE','diagram').replace(',', ' ')}"])
    w.writerow([f"# column x: {meta.get('XTEXT','x')}   column y: {meta.get('YTEXT','y')}"])
    w.writerow(["# boundary = number of the continuous line"])
    w.writerow(["# phases: E = appears or disappears / F0 = fixed at zero amount"])
    w.writerow(["# CALCULATED = boundary obtained by common tangent (not exported from the MAP)"])
    w.writerow(["boundary", "phases", meta.get("XTEXT", "x"), meta.get("YTEXT", "y")])
    w.writerows(filas)

print(f"{ENT} -> {SAL}")
print(f"  {bloque} boundaries, {len(filas)} points")
print(f"  axes: {meta.get('XTEXT')} / {meta.get('YTEXT')}")
