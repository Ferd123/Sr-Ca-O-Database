import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

# Interaction parameters from CaSrO_opt.tdb
def L0(T): return 23756.0123 - 3.61298158 * T
def L1(T): return 916.243186

R = 8.31451

# Liquidus / Solidus parameters (ideal/near-ideal liquidus between CaO: 3222 K, SrO: 2870 K)
T_melt_CaO = 3222.0
T_melt_SrO = 2870.0
dH_melt_CaO = 79500.0
dH_melt_SrO = 70000.0

def solve_misc_gap(T):
    l0 = L0(T)
    l1 = L1(T)
    def dG(y):
        return R * T * np.log(y / (1 - y)) + (1 - 2 * y) * l0 + l1 * (1 - 6 * y + 6 * y**2)
    def mu_CaO(y):
        G = R * T * ((1-y)*np.log(1-y) + y*np.log(y)) + y*(1-y)*(l0 + l1*(1-2*y))
        return G - y * dG(y)
    def mu_SrO(y):
        G = R * T * ((1-y)*np.log(1-y) + y*np.log(y)) + y*(1-y)*(l0 + l1*(1-2*y))
        return G + (1 - y) * dG(y)

    def sys(p):
        return [mu_CaO(p[0]) - mu_CaO(p[1]), mu_SrO(p[0]) - mu_SrO(p[1])]

    sol = root(sys, [0.05, 0.95])
    if sol.success and 0 < sol.x[0] < 0.5 and 0.5 < sol.x[1] < 1:
        return sol.x[0], sol.x[1]
    return None, None

# Compute miscibility gap from 300 K to 1170.4 K
temps = np.linspace(300, 1170.4, 400)
x_cao_b1 = []
x_cao_b2 = []
valid_t = []

for T in temps:
    y1, y2 = solve_misc_gap(T)
    if y1 is not None and abs(y1 - y2) > 0.0005:
        x_cao_b1.append(1 - y1)
        x_cao_b2.append(1 - y2)
        valid_t.append(T)

# Approximate liquidus & solidus curves
x_l = np.linspace(0, 1, 200)
t_liquidus = T_melt_SrO * (1 - x_l) + T_melt_CaO * x_l
t_solidus = t_liquidus - 12 * np.sin(np.pi * x_l)

# Experimental points (Risold 1997 / Jacob 2000)
exp_data = [
    (880, 0.114), (918, 0.906), (937, 0.896), (957, 0.885), (1006, 0.203),
    (1055, 0.817), (1067, 0.805), (1081, 0.790), (1097, 0.772), (1113, 0.737),
    (1117, 0.336), (1125, 0.355), (1128, 0.727), (1100, 0.760), (1100, 0.288)
]
exp_t = [pt[0] for pt in exp_data]
exp_x = [pt[1] for pt in exp_data]

plt.figure(figsize=(10, 7), dpi=300)

# Liquidus & Solidus
plt.plot(x_l, t_liquidus, 'r-', linewidth=2, label='Líquido / Líquidus')
plt.plot(x_l, t_solidus, 'r--', linewidth=1.8, label='Sólidus Halita')

# Miscibility gap
plt.plot(x_cao_b1, valid_t, 'b-', linewidth=2.2, label='Binodal Halita (300 K - 1170 K)')
plt.plot(x_cao_b2, valid_t, 'b-', linewidth=2.2)

# Region text labels
plt.text(0.5, 2000, 'Solución Sólida Homogénea\nHalita (Ca,Sr)O', fontsize=12, ha='center', bbox=dict(boxstyle='round,pad=0.5', facecolor='white', alpha=0.8))
plt.text(0.5, 700, 'Laguna de Miscibilidad\nHalita#1 + Halita#2', fontsize=11, ha='center', color='darkblue', bbox=dict(boxstyle='round,pad=0.4', facecolor='aliceblue', alpha=0.9))
plt.text(0.5, 3150, 'LÍQUIDO', fontsize=12, fontweight='bold', ha='center', color='darkred')

# Experimental points
plt.scatter(exp_x, exp_t, color='crimson', s=45, zorder=5, label='Datos Experimentales (880-1128 K)')

plt.title('Diagrama de Fases Completo CaO-SrO (300 K a 3400 K)\nBase Optimizada CaSrO_opt.tdb', fontsize=14, fontweight='bold', pad=12)
plt.xlabel('Fracción Molar de CaO, $x_{\mathrm{CaO}}$', fontsize=12)
plt.ylabel('Temperatura $T$ (K)', fontsize=12)

plt.xlim(0, 1)
plt.ylim(300, 3400)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10, loc='center left')

plt.tight_layout()
plt.savefig('diagrama_CaSrO_completo.png')
print('Full diagram saved as diagrama_CaSrO_completo.png')
