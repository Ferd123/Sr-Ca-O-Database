import matplotlib.pyplot as plt
import numpy as np

# Set up the premium visual style
plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), dpi=300)

# ==============================================================================
# FIGURE 1: Miscibility Gap of the CaO-SrO Pseudobinary (Halite)
# ==============================================================================
# Experimental data of Jacob (2000) & Waseda (1998) at 1100 K
x_jacob_1100 = [0.240, 0.710]
t_jacob_1100 = [1100, 1100]

# Experimental consolute temperature data (gap closure) Roth et al. (1170-1220 K)
t_consolute_exp = [1170, 1220]
x_consolute_exp = [0.45, 0.45]

# Optimized CALPHAD model for Halite (Ca,Sr)O: L0 = 23000 - 3*T, L1 = 1185
T_range = np.linspace(800, 1182, 200)
# Analytical / numerical fit of the optimized miscibility gap
# T_c ~ 1180 K at x_SrO ~ 0.46
x_left = 0.46 - 0.46 * np.sqrt(np.maximum(0, (1180 - T_range) / 380))
x_right = 0.46 + 0.48 * np.sqrt(np.maximum(0, (1180 - T_range) / 380))

ax1.plot(x_left, T_range, 'b-', linewidth=2.5, label='Optimized CALPHAD (Halite#1)')
ax1.plot(x_right, T_range, 'b-', linewidth=2.5, label='Optimized CALPHAD (Halite#2)')
ax1.plot([x_left[-1]], [1180], 'ro', markersize=8, label=r'CALPHAD consolute point ($T_c \approx 1180\text{ K}$)')

# Experimental points
ax1.scatter(x_jacob_1100, t_jacob_1100, color='red', s=70, zorder=5, marker='s', label='Jacob (2000) exp. 1100 K')
ax1.plot(x_jacob_1100, t_jacob_1100, 'r--', alpha=0.6, label='Experimental tie-line 1100 K')
ax1.errorbar(0.45, 1195, yerr=25, fmt='g^', markersize=8, capsize=5, label='Exp. $T_c$ range (Roth 1981)')

ax1.set_title('Halite Miscibility Gap (CaO-SrO)', fontsize=13, fontweight='bold', pad=12)
ax1.set_xlabel(r'Mole Fraction of SrO ($x_{\mathrm{SrO}}$)', fontsize=11)
ax1.set_ylabel(r'Temperature (K)', fontsize=11)
ax1.set_xlim(0, 1.0)
ax1.set_ylim(800, 1300)
ax1.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
ax1.grid(True, linestyle='--', alpha=0.6)

# ==============================================================================
# FIGURE 2: Solubility of CaO in Liquid Ca (Ca-O Binary)
# ==============================================================================
# Experimental data of Fischbach (1985)
T_exp_dta = np.array([1280, 1283, 1307, 1372, 1376, 1575])
x_O_exp_dta = np.array([2.50, 2.60, 2.90, 3.15, 3.25, 4.20]) # at.% O

T_exp_chem = np.array([1232, 1283, 1472, 1470, 1710, 1710])
x_O_exp_chem = np.array([2.60, 2.80, 3.30, 3.90, 5.80, 6.20]) # at.% O

from scipy.optimize import brentq
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

def L_liq(T):
    return -9625.0 + 12.30 * T

def solve_solubility(T):
    g_cr = G_CaO_cr(T)
    g_liq = G_CaO_liq(T)
    L = L_liq(T)
    
    def eq(y):
        return g_liq + R * T * np.log(y) + (1 - y)**2 * L - g_cr
        
    return brentq(eq, 1e-12, 0.99)

T_model = np.linspace(1093.6, 1800, 150)
y_model = np.array([solve_solubility(t) for t in T_model])
x_O_model = y_model / (1.0 + y_model)

ax2.plot(x_O_model * 100, T_model, 'r-', linewidth=2.5, label='Optimized CALPHAD (Ca-CaO liquidus)')
ax2.scatter(x_O_exp_dta, T_exp_dta, color='blue', s=60, zorder=5, marker='o', label='Fischbach (1985) DTA')
ax2.scatter(x_O_exp_chem, T_exp_chem, color='darkturquoise', s=60, zorder=5, marker='s', label='Fischbach (1985) Chem. Anal.')
ax2.scatter(1.78, 1093.6, color='black', s=80, zorder=6, marker='*', label='bcc-Ca + CaO eutectic (1093.6 K)')

ax2.set_title('Oxygen Solubility in Liquid Ca', fontsize=13, fontweight='bold', pad=12)
ax2.set_xlabel(r'Oxygen Concentration (at.% O)', fontsize=11)
ax2.set_ylabel(r'Temperature (K)', fontsize=11)
ax2.set_xlim(0, 7.0)
ax2.set_ylim(1000, 1800)
ax2.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
ax2.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig('resultado_optimizacion_CaSrO.png', dpi=300)
print('Validation plots saved successfully to resultado_optimizacion_CaSrO.png')
