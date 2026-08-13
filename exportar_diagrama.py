r"""
Convierte un archivo DATAPLOT (.exp) exportado por Thermo-Calc a un CSV
ordenado, apto para ser graficado con cualquier herramienta.

    python exportar_diagrama.py dos.exp diagrama_CaSrO.csv

El .exp contiene el diagrama ya calculado como coordenadas numericas,
agrupadas en bloques. Cada bloque es una linea continua del diagrama y va
precedido de comentarios que declaran las fases implicadas:

    $ BLOCK #68 1 FOR:
    $E LIQUID          <- fase que aparece o desaparece
    $F0 HALITE#2       <- fase fijada en cantidad cero (el limite)

El CSV resultante lleva una fila por punto y una columna que identifica a
que frontera pertenece, de modo que quien lo reciba pueda trazar cada linea
por separado sin disponer de Thermo-Calc.
"""
import csv
import io
import re
import sys

ENT = sys.argv[1] if len(sys.argv) > 1 else "dos.exp"
SAL = sys.argv[2] if len(sys.argv) > 2 else "diagrama_CaSrO.csv"

lineas = io.open(ENT, encoding="latin-1").read().splitlines()

meta = {}
for l in lineas:                       # encabezado global
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
    m = re.match(r"\$([EF]\d*)\s+(\S+)", s)          # $E / $F0 declaran fases
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
    # Las lineas de comentario no deben contener comas: el modulo csv las
    # entrecomillaria y dejarian de reconocerse como comentario al leer el
    # archivo con pandas u hojas de calculo.
    w.writerow([f"# {meta.get('TITLE','diagrama').replace(',', ' ')}"])
    w.writerow([f"# columna x: {meta.get('XTEXT','x')}   columna y: {meta.get('YTEXT','y')}"])
    w.writerow(["# frontera = numero de linea continua"])
    w.writerow(["# fases: E = aparece o desaparece / F0 = fijada en cantidad cero"])
    w.writerow(["# CALCULADO = frontera obtenida por tangente comun (no exportada del MAP)"])
    w.writerow(["frontera", "fases", meta.get("XTEXT", "x"), meta.get("YTEXT", "y")])
    w.writerows(filas)

print(f"{ENT} -> {SAL}")
print(f"  {bloque} fronteras, {len(filas)} puntos")
print(f"  ejes: {meta.get('XTEXT')} / {meta.get('YTEXT')}")
