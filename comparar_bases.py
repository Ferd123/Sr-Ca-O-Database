r"""
Diagrama pseudobinario CaO-SrO calculado con las dos bases, para ver el
efecto del unico parametro que las diferencia: L(LIQUID,CAO,SRO;0).

Izquierda  CaSrO_opt.tdb      L = +25000 (Risold 1997)
Derecha    CaSrO_opt_liq.tdb  L =  +2000 (analogia Zhang 2016)

La laguna de la halita es identica en ambas: no depende de ese parametro.
"""
import numpy as np
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt

exec(open("liquidus_CaO_SrO.py", encoding="utf-8").read().split("L0_nuestro =")[0])

R = 8.31451
L0 = lambda T: 23756.0 - 3.6130 * T
L1 = 916.24


def binodal_halita(T, guess=(0.05, 0.95)):
    g = lambda y: G_sol(y, T, L0, L1)
    d = lambda y: (g(y + 1e-7) - g(y - 1e-7)) / 2e-7

    def eqs(p):
        y1, y2 = p
        if not (1e-9 < y1 < 1 - 1e-9 and 1e-9 < y2 < 1 - 1e-9):
            return [1e3, 1e3]
        m = (g(y2) - g(y1)) / (y2 - y1)
        return [d(y1) - m, d(y2) - m]

    p, _, ier, _ = fsolve(eqs, guess, full_output=True)
    return None if (ier != 1 or abs(p[0] - p[1]) < 1e-4) else sorted(p)


def critico():
    def T_de_y(y):
        d = 1 / (1 - y) ** 2 - 1 / y ** 2
        return -12 * L1 / (R * d) if abs(d) > 1e-12 else None

    def f(y):
        T = T_de_y(y)
        if T is None or T <= 0:
            return 1e9
        return R * T * (1 / (1 - y) + 1 / y) - 2 * (L0(T) + 3 * L1) + 12 * L1 * y

    y = brentq(f, 0.35, 0.499999)
    return y, T_de_y(y)


YC, TC = critico()
Ts_liq = np.linspace(2840, 3221.5, 600)

fig, axes = plt.subplots(1, 2, figsize=(12.4, 6.2), sharey=True)
for ax, (L, tit, col) in zip(axes, [(25000.0, "CaSrO_opt.tdb\nL(liq) = +25000   (Risold 1997)", "#2a78d6"),
                                    (2000.0, "CaSrO_opt_liq.tdb\nL(liq) = +2000   (analogia Zhang 2016)", "#1baf7a")]):
    T, ys, yl = campo(L, L0, L1, Ts_liq)
    ax.plot(yl, T, ".", color=col, ms=2.4)
    ax.plot(ys, T, ".", color=col, ms=1.2, alpha=0.6)

    # laguna de la halita, identica en ambas
    r1, r2, Tg = [], [], []
    g = (0.02, 0.98)
    for t in np.linspace(400, TC, 400)[:-1]:
        b = binodal_halita(t, g)
        if b is None:
            continue
        g = tuple(b); Tg.append(t); r1.append(b[0]); r2.append(b[1])
    ax.plot(r1 + [YC] + r2[::-1], Tg + [TC] + Tg[::-1], color="#eb6834", lw=1.9)
    ax.plot([YC], [TC], "*", ms=11, color="#0b0b0b", mec="white", mew=0.8)

    ax.set_title(tit, fontsize=10.5)
    ax.set_xlabel("fraccion de SrO     (0 = CaO,  1 = SrO)")
    ax.grid(color="#e1e0d9", lw=0.7); ax.set_axisbelow(True)
    ax.set_xlim(0, 1); ax.set_ylim(400, 3300)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)

axes[0].set_ylabel("T / K")
axes[0].plot([], [], color="#eb6834", lw=1.9, label=f"laguna de la halita (Tc = {TC:.0f} K, identica)")
axes[0].plot([], [], ".", color="#52514e", ms=6, label="liquidus / solidus")
axes[0].legend(frameon=False, fontsize=8.5, loc="center left")
fig.suptitle("Un solo parametro de diferencia entre las dos bases", fontsize=12, y=0.98)
fig.tight_layout()
fig.savefig("comparacion_bases.png", dpi=200)
print(f"Tc halita = {TC:.2f} K,  x(CaO)c = {1-YC:.4f}   (igual en ambas)")
print("-> comparacion_bases.png")
