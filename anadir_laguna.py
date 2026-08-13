r"""
Appends to the diagram CSV the boundary of the halite miscibility gap, which
the Thermo-Calc MAP command does not trace automatically.

    python anadir_laguna.py diagrama_CaSrO_completo.csv

The binodal is obtained by solving the common tangent condition on the model of
the optimized database:

    G = RT[(1-y)ln(1-y) + y ln y] + (1-y)y[L0 + L1(1-2y)]
    L0 = V1 + V2*T,  L1 = V3
    V1 = 23756.012326,  V2 = -3.6129815828,  V3 = 916.2431862683

y = site fraction of Sr on the cation sublattice. The CSV axis is the
Thermo-Calc mole fraction, X(SR) = y/2, because it counts the oxygen as well.

The appended rows are tagged CALCULATED so that it is on record that they do
not come from the MAP export.
"""
import csv
import io
import sys

import numpy as np
from scipy.optimize import fsolve

CSV = sys.argv[1] if len(sys.argv) > 1 else "diagrama_CaSrO_completo.csv"
R = 8.31451
V1, V2, V3 = 23756.012326, -3.6129815828, 916.2431862683


def g(y, T):
    x = 1.0 - y
    return R * T * (x * np.log(x) + y * np.log(y)) + x * y * ((V1 + V2 * T) + V3 * (x - y))


def binodal(T, guess=(0.10, 0.90)):
    d = lambda y: (g(y + 1e-7, T) - g(y - 1e-7, T)) / 2e-7

    def eqs(p):
        y1, y2 = p
        if not (1e-9 < y1 < 1 - 1e-9 and 1e-9 < y2 < 1 - 1e-9):
            return [1e3, 1e3]
        m = (g(y2, T) - g(y1, T)) / (y2 - y1)
        return [d(y1) - m, d(y2) - m]

    p, _, ier, _ = fsolve(eqs, guess, full_output=True)
    if ier != 1 or abs(p[0] - p[1]) < 1e-4:
        return None
    return sorted(p)


def critico():
    """Exact critical point: d2G/dy2 = d3G/dy3 = 0.

    With L0 = V1 + V2*T and L1 = V3 the derivatives are analytic:
        d2G = RT[1/(1-y) + 1/y] - 2(L0 + 3*L1) + 12*L1*y
        d3G = RT[1/(1-y)^2 - 1/y^2] + 12*L1
    T is solved from d3G = 0 and substituted into d2G = 0, which then has a
    single unknown. Bisecting on the failure of fsolve underestimates the
    critical temperature by a little over 1 K, because fsolve stops converging
    somewhat before the top.
    """
    from scipy.optimize import brentq

    def T_de_y(y):
        den = 1 / (1 - y) ** 2 - 1 / y ** 2
        return -12 * V3 / (R * den) if abs(den) > 1e-12 else None

    def f(y):
        T = T_de_y(y)
        if T is None or T <= 0:
            return 1e9
        return R * T * (1 / (1 - y) + 1 / y) - 2 * ((V1 + V2 * T) + 3 * V3) + 12 * V3 * y

    y = brentq(f, 0.35, 0.499999)
    return y, T_de_y(y)


Y_C, TC = critico()

filas = list(csv.reader(io.open(CSV, encoding="utf-8")))
nueva = max(int(r[0]) for r in filas if r and r[0].isdigit()) + 1

# Continuation: each solution serves as the initial guess for the next, which
# allows getting much closer to the top than always starting from (0.1, 0.9).
rama1, rama2 = [], []
guess = (0.02, 0.98)
for T in np.linspace(300.0, TC, 500)[:-1]:
    b = binodal(T, guess)
    if b is None:
        continue
    guess = tuple(b)
    rama1.append((b[0] / 2, T))
    rama2.append((b[1] / 2, T))

# Exact closure at the analytically computed critical point.
curva = rama1 + [(Y_C / 2, TC)] + list(reversed(rama2))
etiqueta = "CALCULATED E:HALITE#1 + E:HALITE#2 (miscibility gap)"

with io.open(CSV, "a", encoding="utf-8", newline="") as f:
    w = csv.writer(f)
    for x, T in curva:
        w.writerow([nueva, etiqueta, f"{x:.10f}", f"{T:.4f}"])

print(f"appended boundary {nueva} to file {CSV}")
print(f"  {len(curva)} points")
print(f"  critical temperature {TC:.2f} K = {TC - 273.15:.2f} C")
print(f"  critical X(SR) {Y_C / 2:.4f}")
