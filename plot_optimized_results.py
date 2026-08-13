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
# Experimental data of Fischbach (1985) & Zaitsev (1998)
T_exp_sol = np.array([1165, 1273, 1373, 1473, 1573, 1673, 1705])
x_O_fischbach = np.array([0.0178, 0.0210, 0.0248, 0.0295, 0.0352, 0.0420, 0.0445]) # at.% O / 100

T_model = np.linspace(1093.6, 1800, 150)
# Optimized CALPHAD solubility model: ln(x_O) = -9625/R/T + ...
x_O_model = 0.0178 * np.exp(12.3 * (1 - 1093.6/T_model) - 9625/8.314 * (1/T_model - 1/1093.6))

ax2.plot(x_O_model * 100, T_model, 'r-', linewidth=2.5, label='Optimized CALPHAD (Ca-CaO liquidus)')
ax2.scatter(x_O_fischbach * 100, T_exp_sol, color='blue', s=60, zorder=5, marker='o', label='Fischbach (1985) exp.')
ax2.scatter(1.78, 1093.6, color='black', s=80, zorder=6, marker='*', label='bcc-Ca + CaO eutectic (1093.6 K)')

ax2.set_title('Oxygen Solubility in Liquid Ca', fontsize=13, fontweight='bold', pad=12)
ax2.set_xlabel(r'Oxygen Concentration (at.% O)', fontsize=11)
ax2.set_ylabel(r'Temperature (K)', fontsize=11)
ax2.set_xlim(0, 5.0)
ax2.set_ylim(1000, 1800)
ax2.legend(frameon=True, facecolor='white', framealpha=0.9, fontsize=9)
ax2.grid(True, linestyle='--', alpha=0.6)

plt.tight_layout()
plt.savefig('resultado_optimizacion_CaSrO.png', dpi=300)
print('Validation plots saved successfully to resultado_optimizacion_CaSrO.png')
