import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import root

# Interaction parameters from CaSrO_opt.tdb
# L0 = 23756.01 - 3.61298 * T
# L1 = 916.24
R = 8.31451

def L0(T): return 23756.0123 - 3.61298158 * T
def L1(T): return 916.243186

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

    sol = root(sys, [0.1, 0.9])
    if sol.success and 0 < sol.x[0] < 0.5 and 0.5 < sol.x[1] < 1:
        return sol.x[0], sol.x[1]
    return None, None

# Compute miscibility gap curve
temps = np.linspace(800, 1170.4, 300)
x_cao_branch1 = []
x_cao_branch2 = []
valid_t = []

for T in temps:
    y1, y2 = solve_misc_gap(T)
    if y1 is not None and abs(y1 - y2) > 0.0005:
        # Convert site fraction of Sr on halite (y_Sr) to mole fraction of CaO: x(CaO) = 1 - y_Sr
        x_cao_branch1.append(1 - y1)
        x_cao_branch2.append(1 - y2)
        valid_t.append(T)

# Experimental data points from Risold et al. (1997) & Jacob et al. (2000)
exp_data_AGAPT = [
    (880, 0.114), (918, 0.906), (937, 0.896), (957, 0.885), (1006, 0.203),
    (1055, 0.817), (1067, 0.805), (1081, 0.790), (1097, 0.772), (1113, 0.737),
    (1117, 0.336), (1125, 0.355), (1128, 0.727)
]
exp_t = [pt[0] for pt in exp_data_AGAPT]
exp_x_cao = [pt[1] for pt in exp_data_AGAPT]

# Experimental point AGAP (1100 K)
exp_agap_t = [1100, 1100]
exp_agap_x_cao = [1 - 0.240, 1 - 0.712] # x(CaO) = 0.760 y 0.288

# Plotting
plt.figure(figsize=(9, 6.5), dpi=300)

# Miscibility gap curve
plt.plot(x_cao_branch1, valid_t, 'b-', linewidth=2.2, label='Calculated binodal (CaSrO_opt.tdb)')
plt.plot(x_cao_branch2, valid_t, 'b-', linewidth=2.2)

# Experimental points
plt.scatter(exp_x_cao, exp_t, color='red', s=45, zorder=5, label='AGAPT data (Risold 1997)')
plt.scatter(exp_agap_x_cao, exp_agap_t, color='darkgreen', marker='s', s=55, zorder=5, label='AGAP data (Jacob 2000)')

plt.title('CaO-SrO Pseudobinary Diagram\nMiscibility Gap of the Halite (Ca,Sr)O', fontsize=14, fontweight='bold', pad=12)
plt.xlabel(r'Mole Fraction of CaO, $x_{\mathrm{CaO}}$', fontsize=12)
plt.ylabel('Temperature $T$ (K)', fontsize=12)

plt.xlim(0, 1)
plt.ylim(800, 1250)
plt.grid(True, linestyle='--', alpha=0.5)
plt.legend(fontsize=10, loc='upper center')

plt.tight_layout()
plt.savefig('diagrama_CaSrO_optimizado.png')
print('Diagram saved as diagrama_CaSrO_optimizado.png')
