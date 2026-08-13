r"""
Actividad del SrO en el fundido CaO-SrO a 3000 K, para los tres valores del
parametro de interaccion discutidos. Muestra que el parametro de exceso ES el
coeficiente de actividad a dilucion infinita, es decir la constante de Henry.

Liquido regular de un solo parametro:
    RT ln(gamma_SrO) = L (1-y)^2      ->  gamma_inf = exp(L/RT)
    a_SrO = y exp[ L (1-y)^2 / RT ]
"""
import numpy as np
import matplotlib.pyplot as plt

R, T = 8.31451, 3000.0
y = np.linspace(1e-4, 1, 400)

casos = [(25000.0, "L = +25000  (base actual)", "#2a78d6"),
         (2000.0,  "L = +2000   (analogia Zhang)", "#eda100"),
         (0.0,     "L = 0       (ideal, Raoult puro)", "#eb6834")]

fig, ax = plt.subplots(figsize=(7.6, 6.2))
ax.plot(y, y, color="#52514e", lw=1.3, ls=(0, (6, 3)), label="ley de Raoult  (a = y)")

for L, et, c in casos:
    a = y * np.exp(L * (1 - y) ** 2 / (R * T))
    ginf = np.exp(L / (R * T))
    ax.plot(y, a, color=c, lw=2.1, label=et)
    # tangente de Henry en el origen: pendiente = gamma_inf
    yy = np.linspace(0, 0.42, 2)
    ax.plot(yy, ginf * yy, color=c, lw=1.0, ls=":", alpha=0.9)
    ax.annotate(f"$\gamma^\infty$ = {ginf:.2f}", (0.42, ginf * 0.42),
                color=c, fontsize=9, va="center", ha="left")

ax.set_xlabel("$y_{SrO}$  en el fundido")
ax.set_ylabel("actividad del SrO")
ax.set_title("El parametro de exceso es la constante de Henry   ·   3000 K")
ax.set_xlim(0, 1.15); ax.set_ylim(0, 1.55)
ax.grid(color="#e1e0d9", lw=0.7); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.text(0.60, 0.10, "punteado = recta de Henry\n(pendiente $\gamma^\infty$)",
        fontsize=8.5, color="#52514e")
fig.tight_layout(); fig.savefig("actividad_raoult_henry.png", dpi=200)
print("-> actividad_raoult_henry.png")
