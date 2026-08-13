r"""
Activity of SrO in the CaO-SrO melt at 3000 K, for the three values of the
interaction parameter under discussion. It shows that the excess parameter IS
the activity coefficient at infinite dilution, that is, the Henry constant.

Regular liquid with a single parameter:
    RT ln(gamma_SrO) = L (1-y)^2      ->  gamma_inf = exp(L/RT)
    a_SrO = y exp[ L (1-y)^2 / RT ]
"""
import numpy as np
import matplotlib.pyplot as plt

R, T = 8.31451, 3000.0
y = np.linspace(1e-4, 1, 400)

casos = [(25000.0, "L = +25000  (current database)", "#2a78d6"),
         (2000.0,  "L = +2000   (Zhang analogy)", "#eda100"),
         (0.0,     "L = 0       (ideal, pure Raoult)", "#eb6834")]

fig, ax = plt.subplots(figsize=(7.6, 6.2))
ax.plot(y, y, color="#52514e", lw=1.3, ls=(0, (6, 3)), label="Raoult's law  (a = y)")

for L, et, c in casos:
    a = y * np.exp(L * (1 - y) ** 2 / (R * T))
    ginf = np.exp(L / (R * T))
    ax.plot(y, a, color=c, lw=2.1, label=et)
    # Henry tangent at the origin: slope = gamma_inf
    yy = np.linspace(0, 0.42, 2)
    ax.plot(yy, ginf * yy, color=c, lw=1.0, ls=":", alpha=0.9)
    ax.annotate(rf"$\gamma^\infty$ = {ginf:.2f}", (0.42, ginf * 0.42),
                color=c, fontsize=9, va="center", ha="left")

ax.set_xlabel("$y_{SrO}$  in the melt")
ax.set_ylabel("activity of SrO")
ax.set_title("The excess parameter is the Henry constant   ·   3000 K")
ax.set_xlim(0, 1.15); ax.set_ylim(0, 1.55)
ax.grid(color="#e1e0d9", lw=0.7); ax.set_axisbelow(True)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=9, loc="upper left")
ax.text(0.60, 0.10, "dotted = Henry line\n" r"(slope $\gamma^\infty$)",
        fontsize=8.5, color="#52514e")
fig.tight_layout(); fig.savefig("actividad_raoult_henry.png", dpi=200)
print("-> actividad_raoult_henry.png")
