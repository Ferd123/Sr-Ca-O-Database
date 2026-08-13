"""
Energia de Gibbs de mezcla de la halita (Ca,Sr)O y construccion de la
tangente comun que define la laguna de miscibilidad.

Modelo (CaSrO.TDB, parametros de Risold et al. 1997, Tabla VI):

    dGmix = R*T*(y_Ca*ln y_Ca + y_Sr*ln y_Sr)
          + y_Ca*y_Sr*(L0 + L1*(y_Ca - y_Sr))
    L0 = 23000 - 3*T        L1 = 1185      [J/mol de (Ca,Sr)O]

y = fraccion de sitio en la primera subred. La fraccion molar de
Thermo-Calc es X(SR) = y_Sr/2, porque cuenta tambien el oxigeno.

Genera dG_halita_laguna.png
"""

import numpy as np
from scipy.optimize import fsolve, brentq
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

R = 8.31451

def L0(T): return 23000.0 - 3.0 * T
def L1(T): return 1185.0

def dgmix(y, T):
    """Energia de Gibbs de mezcla, J por mol de (Ca,Sr)O. y = y_Sr."""
    x = 1.0 - y
    return R * T * (x * np.log(x) + y * np.log(y)) + x * y * (L0(T) + L1(T) * (x - y))

def d1(y, T, h=1e-7):
    return (dgmix(y + h, T) - dgmix(y - h, T)) / (2 * h)

def d2(y, T, h=1e-5):
    return (dgmix(y + h, T) - 2 * dgmix(y, T) + dgmix(y - h, T)) / h ** 2

def binodal(T):
    """Tangente comun: dG/dy igual en los dos puntos e igual a la secante."""
    def eqs(p):
        y1, y2 = p
        if not (1e-9 < y1 < 1 - 1e-9 and 1e-9 < y2 < 1 - 1e-9):
            return [1e3, 1e3]
        m = (dgmix(y2, T) - dgmix(y1, T)) / (y2 - y1)
        return [d1(y1, T) - m, d1(y2, T) - m]
    p, _, ier, _ = fsolve(eqs, [0.10, 0.90], full_output=True)
    if ier != 1 or abs(p[0] - p[1]) < 1e-4:
        return None
    return tuple(sorted(p))

def spinodal(T):
    """Puntos de inflexion, d2G/dy2 = 0."""
    ys = np.linspace(0.001, 0.999, 4000)
    v = d2(ys, T)
    cruces = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    if len(cruces) < 2:
        return None
    return (brentq(d2, ys[cruces[0]], ys[cruces[0] + 1], args=(T,)),
            brentq(d2, ys[cruces[-1]], ys[cruces[-1] + 1], args=(T,)))

# temperatura critica
lo, hi = 1000.0, 2000.0
for _ in range(200):
    mid = (lo + hi) / 2
    if binodal(mid) is None:
        hi = mid
    else:
        lo = mid
TC = lo

# datos experimentales, X(SR) -> y_Sr = 2*X(SR)
EXP = [(0.4430, 880), (0.0470, 918), (0.0520, 937), (0.0575, 957),
       (0.3985, 1006), (0.0915, 1055), (0.0975, 1067), (0.1050, 1081),
       (0.1140, 1097), (0.1315, 1113), (0.3320, 1117), (0.3225, 1125),
       (0.1365, 1128)]
EXP_1100 = [(0.120, 1100), (0.356, 1100)]

# ---------------------------------------------------------------- colores
SURF, INK, INK2, MUTED = "#fcfcfb", "#0b0b0b", "#52514e", "#898781"
GRID, AXIS = "#e1e0d9", "#c3c2b7"
AZUL, NARANJA = "#2a78d6", "#eb6834"
RAMPA = {880: "#86b6ef", 1000: "#5598e7", 1100: "#2a78d6", 1179: "#0d366b"}
TS = [880, 1000, 1100, 1179]
DESTACADA = 1100

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Segoe UI", "DejaVu Sans"],
    "font.size": 9, "axes.linewidth": 0.8,
    "figure.facecolor": SURF, "axes.facecolor": SURF,
})

fig, (ax1, ax2) = plt.subplots(
    2, 1, figsize=(7.2, 8.4), sharex=True,
    gridspec_kw={"height_ratios": [1.15, 1], "hspace": 0.10})

# ================================================================ panel A
y = np.linspace(1e-4, 1 - 1e-4, 1200)

for T in TS:
    ax1.plot(y, dgmix(y, T) / 1000, color=RAMPA[T],
             lw=2.4 if T == DESTACADA else 1.6, zorder=3, solid_capstyle="round")
    b = binodal(T)
    if b is None:
        continue
    y1, y2 = b
    g1, g2 = dgmix(y1, T) / 1000, dgmix(y2, T) / 1000
    ax1.plot([y1, y2], [g1, g2], color=MUTED, lw=1.4 if T == DESTACADA else 1.0,
             ls=(0, (5, 3)), zorder=2)
    ax1.plot([y1, y2], [g1, g2], "o", ms=7 if T == DESTACADA else 5.5,
             color=INK, mec=SURF, mew=1.6, zorder=6)
    s = spinodal(T)
    if s:
        ax1.plot(list(s), [dgmix(v, T) / 1000 for v in s], "o", ms=6,
                 mfc=SURF, mec=NARANJA, mew=1.8, zorder=5)

# etiquetas directas de temperatura, cada una sobre su propia curva
for T, xa, dy, va in [(880, 0.46, 0.06, "bottom"), (1000, 0.46, -0.04, "top"),
                      (1100, 0.46, 0.07, "bottom"), (1179, 0.50, -0.10, "top")]:
    ax1.annotate(f"{T} K" + "  (Tc)" * (T == 1179),
                 (xa, dgmix(xa, T) / 1000 + dy), color=RAMPA[T], va=va,
                 fontsize=8.5, fontweight="bold", ha="center", zorder=7)

# nota explicativa, en la esquina que queda libre
ax1.text(0.028, 1.95,
         "La tangente comun toca la curva en dos puntos.\n"
         "Esas dos composiciones tienen el mismo potencial\n"
         "quimico, asi que coexisten: ahi esta la laguna.\n"
         "Al subir T la curva se aplana, los dos puntos se\n"
         "acercan, y en Tc se funden en uno solo.",
         fontsize=8.2, color=INK2, va="top", linespacing=1.55, zorder=7)

ax1.axhline(0, color=AXIS, lw=0.8, zorder=1)
ax1.set_ylabel("$\\Delta G_{mezcla}$   /  kJ por mol de (Ca,Sr)O", color=INK2)
ax1.set_ylim(-2.6, 2.15)
ax1.set_title("Halita (Ca,Sr)O:  la tangente comun y la laguna que produce",
              loc="left", fontsize=12, fontweight="bold", color=INK, pad=24)
ax1.text(0, 1.022, "Parametros de Risold et al. 1997 (CaSrO.TDB).  "
         "Puntos llenos: binodal.  Huecos: espinodal.",
         transform=ax1.transAxes, fontsize=8.5, color=MUTED, va="bottom")

leg = [Line2D([], [], color=AZUL, lw=2.2, label="$\\Delta G_{mezcla}$"),
       Line2D([], [], color=MUTED, lw=1.3, ls=(0, (5, 3)), label="tangente comun"),
       Line2D([], [], marker="o", ls="", ms=7, color=INK, mec=SURF, mew=1.4,
              label="binodal (composicion de equilibrio)"),
       Line2D([], [], marker="o", ls="", ms=6, mfc=SURF, mec=NARANJA, mew=1.8,
              label="espinodal (limite de estabilidad)")]
ax1.legend(handles=leg, loc="upper right", frameon=False, fontsize=8.2,
           labelcolor=INK2, handletextpad=0.7, borderaxespad=0.4)

# ================================================================ panel B
Tg = np.linspace(300, TC - 0.05, 500)
bs = [binodal(T) for T in Tg]
sp = [spinodal(T) for T in Tg]
ok = [i for i, v in enumerate(bs) if v]
ax2.plot([bs[i][0] for i in ok] + [bs[i][1] for i in reversed(ok)],
         [Tg[i] for i in ok] + [Tg[i] for i in reversed(ok)],
         color=AZUL, lw=2.2, zorder=4, solid_capstyle="round")
oks = [i for i, v in enumerate(sp) if v]
ax2.plot([sp[i][0] for i in oks] + [sp[i][1] for i in reversed(oks)],
         [Tg[i] for i in oks] + [Tg[i] for i in reversed(oks)],
         color=NARANJA, lw=1.6, ls=(0, (5, 3)), zorder=3)

ax2.fill(([bs[i][0] for i in ok] + [bs[i][1] for i in reversed(ok)]),
         ([Tg[i] for i in ok] + [Tg[i] for i in reversed(ok)]),
         color=AZUL, alpha=0.055, lw=0, zorder=1)

ax2.plot([2 * x for x, _ in EXP], [t for _, t in EXP], "s", ms=5.5,
         mfc=SURF, mec=INK, mew=1.3, zorder=6)
ax2.plot([2 * x for x, _ in EXP_1100], [t for _, t in EXP_1100], "o", ms=7,
         color=INK, mec=SURF, mew=1.4, zorder=7)

ax2.plot([0.5], [TC], "*", ms=13, color=INK, mec=SURF, mew=1.0, zorder=8)
ax2.annotate(f"punto critico  {TC:.0f} K", (0.5, TC), xytext=(0.5, TC + 55),
             ha="center", fontsize=8.5, color=INK2)

# lineas guia que conectan los dos paneles
for yy in binodal(DESTACADA):
    ax2.plot([yy, yy], [300, DESTACADA], color=MUTED, lw=0.7, ls=":", zorder=2)
    ax1.plot([yy, yy], [-2.6, dgmix(yy, DESTACADA) / 1000], color=MUTED,
             lw=0.7, ls=":", zorder=0)
ax2.plot(list(binodal(DESTACADA)), [DESTACADA] * 2, "o", ms=7, color=INK,
         mec=SURF, mew=1.4, zorder=9)
ax2.annotate("1100 K", (binodal(DESTACADA)[1] + 0.035, DESTACADA), fontsize=8.5,
             color=INK2, va="center")

ax2.annotate("una sola halita\nhomogenea", (0.5, TC + 210), ha="center",
             fontsize=8.5, color=MUTED)
ax2.annotate("dos halitas:\nuna rica en Ca, otra rica en Sr", (0.5, 640),
             ha="center", fontsize=8.5, color=INK2)

leg2 = [Line2D([], [], color=AZUL, lw=2.2, label="binodal, calculada"),
        Line2D([], [], color=NARANJA, lw=1.6, ls=(0, (5, 3)), label="espinodal, calculada"),
        Line2D([], [], marker="s", ls="", ms=5.5, mfc=SURF, mec=INK, mew=1.3,
               label="Fig. 1 de Risold 1997 (13 pts)"),
        Line2D([], [], marker="o", ls="", ms=7, color=INK, mec=SURF, mew=1.4,
               label="Jacob y Waseda 1998, 1100 K")]
ax2.legend(handles=leg2, loc="lower center", frameon=False, fontsize=8.2,
           labelcolor=INK2, ncol=2, handletextpad=0.7, borderaxespad=0.8)

ax2.set_xlabel("$y_{Sr}$  =  fraccion de sitio de Sr en la subred cationica"
               "        [ X(SR) de Thermo-Calc = $y_{Sr}$/2 ]", color=INK2)
ax2.set_ylabel("Temperatura  /  K", color=INK2)
ax2.set_xlim(0, 1)
ax2.set_ylim(300, TC + 330)

for ax in (ax1, ax2):
    ax.grid(True, color=GRID, lw=0.7, zorder=0)
    ax.set_axisbelow(True)
    for s_ in ("top", "right"):
        ax.spines[s_].set_visible(False)
    for s_ in ("left", "bottom"):
        ax.spines[s_].set_color(AXIS)
    ax.tick_params(colors=MUTED, labelcolor=INK2, length=3)

fig.text(0.008, 0.005, "CaSrO.TDB  ·  L0 = 23000 - 3T,  L1 = 1185 J/mol",
         fontsize=7.5, color=MUTED)
fig.savefig("dG_halita_laguna.png", dpi=200, bbox_inches="tight", facecolor=SURF)
print(f"Tc = {TC:.1f} K")
for T in TS:
    b = binodal(T)
    if b:
        print(f"{T} K: binodal y_Sr = {b[0]:.4f} / {b[1]:.4f}   "
              f"X(SR) = {b[0]/2:.4f} / {b[1]/2:.4f}")
print("-> dG_halita_laguna.png")
