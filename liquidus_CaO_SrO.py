r"""
Liquidus y solidus del pseudobinario CaO-SrO por construccion de envolvente
convexa sobre las curvas de energia de Gibbs del liquido (asociados) y de la
halita (dos subredes), con las funciones de CaSrO_opt.tdb.

Objetivo: mostrar de que depende la forma del campo bifasico L+halita, que es
donde nuestro diagrama difiere del de FactSage FToxid.

La envolvente convexa se usa en lugar de resolver la tangente comun con un
metodo de Newton porque no depende de valores iniciales y detecta cualquier
topologia: lente simple, minimo congruente o eutectico.

Variable de composicion: y = fraccion de SrO. Energias por mol de formula
(Ca,Sr)O, es decir dos atomos.
"""
import numpy as np
import matplotlib.pyplot as plt

R = 8.31451


def gein(theta, T):
    return 1.5 * R * theta + 3 * R * T * np.log(1 - np.exp(-theta / T))


def G_CaO_cr(T):
    T = np.asarray(T, float)
    baja = (-652134.4 + 1.142993 * gein(369.447, T) + 0.62542 * gein(601.229, T)
            + 0.218718 * gein(188.291, T) - 0.00182458 * T ** 2
            - 3.1187482541e-02 * np.exp(0.00304142 * T))
    alta = (-726686.2348 + 556.326025 * T - 79.464890 * T * np.log(T)
            - 2.269506e-04 * T ** 2)
    return np.where(T <= 3222.0, baja, alta)


def G_CaO_liq(T):
    return (-656748.5662 + 572.751716 * T - 84.370735 * T * np.log(T)
            + 2.374069e-04 * T ** 2)


def G_SrO_cr(T):
    return (-607870 + 268.9 * T - 47.56 * T * np.log(T)
            - 0.00307 * T ** 2 + 190000 * T ** (-1))


def G_SrO_liq(T):
    return -566346 + 449.0 * T - 73.1 * T * np.log(T)


def ideal(y, T):
    y = np.clip(y, 1e-12, 1 - 1e-12)
    return R * T * (y * np.log(y) + (1 - y) * np.log(1 - y))


def G_liq(y, T, L_liq):
    return (1 - y) * G_CaO_liq(T) + y * G_SrO_liq(T) + ideal(y, T) + (1 - y) * y * L_liq


def G_sol(y, T, L0, L1):
    return ((1 - y) * G_CaO_cr(T) + y * G_SrO_cr(T) + ideal(y, T)
            + (1 - y) * y * (L0(T) + L1 * ((1 - y) - y)))


def envolvente(T, L_liq, L0, L1, n=1200):
    """Tie-lines a temperatura T por envolvente convexa inferior.

    Devuelve una lista de (y_izq, fase_izq, y_der, fase_der) para cada arista
    de la envolvente que une dos puntos no contiguos, es decir cada region
    bifasica.
    """
    y = np.linspace(1e-5, 1 - 1e-5, n)
    pts = np.concatenate([np.c_[y, G_liq(y, T, L_liq), np.zeros(n)],
                          np.c_[y, G_sol(y, T, L0, L1), np.ones(n)]])
    orden = np.argsort(pts[:, 0], kind="stable")
    pts = pts[orden]

    # Envolvente convexa INFERIOR (cadena monotona de Andrew). Se descarta el
    # ultimo punto mientras el giro no sea antihorario, es decir mientras el
    # producto cruzado sea <= 0.
    def cruz(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    casco = []
    for p in pts:
        while len(casco) >= 2 and cruz(casco[-2], casco[-1], p) <= 0:
            casco.pop()
        casco.append(p)

    tie = []
    for (x1, g1, f1), (x2, g2, f2) in zip(casco[:-1], casco[1:]):
        if x2 - x1 > 3.0 / n and f1 != f2:      # arista larga entre fases distintas
            tie.append((x1, int(f1), x2, int(f2)))
    return tie


def campo(L_liq, L0, L1, Ts):
    """Devuelve (T, y_solido, y_liquido) del campo L+halita."""
    Tl, ysol, yliq = [], [], []
    for T in Ts:
        for x1, f1, x2, f2 in envolvente(T, L_liq, L0, L1):
            ys = x1 if f1 == 1 else x2
            yl = x2 if f1 == 1 else x1
            Tl.append(T); ysol.append(ys); yliq.append(yl)
    return np.array(Tl), np.array(ysol), np.array(yliq)


L0_nuestro = lambda T: 23756.0 - 3.6130 * T
L1_nuestro = 916.24
Ts = np.linspace(2830, 3221.5, 700)

casos = [("L(liq) = +25000   (Risold, la de esta base)", 25000.0, "#2a78d6"),
         ("L(liq) = +10000", 10000.0, "#eda100"),
         ("L(liq) = 0   (liquido ideal)", 0.0, "#eb6834")]

fig, ax = plt.subplots(figsize=(8.2, 6.4))
print(f"{'caso':52s} {'ancho max':>10s} {'en y':>6s} {'T min liquidus':>15s}")
for etiqueta, Ll, color in casos:
    T, ys, yl = campo(Ll, L0_nuestro, L1_nuestro, Ts)
    ax.plot(yl, T, ".", color=color, ms=2.0, label=etiqueta)
    ax.plot(ys, T, ".", color=color, ms=1.0, alpha=0.55)
    anc = np.abs(yl - ys); i = int(np.argmax(anc))
    print(f"{etiqueta:52s} {anc[i]:10.3f} {yl[i]:6.2f} {T.min():13.0f} K")

ax.plot([], [], ".", color="#52514e", ms=6, label="liquidus (trazo grueso)")
ax.plot([], [], ".", color="#52514e", ms=3, alpha=0.55, label="solidus (trazo fino)")
ax.set_xlabel("fraccion de SrO        (0 = CaO,  1 = SrO)")
ax.set_ylabel("T / K")
ax.set_title("El campo L + halita lo fija la DIFERENCIA de exceso entre liquido y solido")
ax.grid(color="#e1e0d9", lw=0.7); ax.set_axisbelow(True)
ax.set_xlim(0, 1); ax.set_ylim(2840, 3245)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.legend(frameon=False, fontsize=8.5, loc="lower left")
fig.tight_layout()
fig.savefig("liquidus_comparacion.png", dpi=200)
print("\n-> liquidus_comparacion.png")
