import numpy as np
import matplotlib.pyplot as plt
import re
import os
import matplotlib as mpl

# =========================================================================
# CONSTANTS & CONFIG (Premium scientific style)
# =========================================================================
mpl.rcParams.update(
    {
        "font.family": "serif",
        "font.serif": ["Times New Roman"],
        "axes.labelsize": 11,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "axes.titlesize": 12,
        "axes.linewidth": 1.0,
        "lines.linewidth": 1.2,
        "xtick.direction": "inout",
        "ytick.direction": "inout",
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "figure.dpi": 600,
        "savefig.dpi": 600,
        "svg.fonttype": "none",
        "mathtext.fontset": "custom",
        "mathtext.rm": "Times New Roman",
        "mathtext.it": "Times New Roman:italic",
        "mathtext.bf": "Times New Roman:bold",
    }
)

# =========================================================================
# 1. PARSE FILE dos.exp
# =========================================================================
exp_path = "dos.exp"
print(f"Reading file: {exp_path}")

try:
    with open(exp_path, "r", encoding="latin-1") as f:
        content = f.read()
except FileNotFoundError:
    print(f"Error: File not found {exp_path}")
    exit()

def parse_thermocalc_exp(content):
    blocks = []
    current_block = None
    number_pattern = re.compile(r"[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?")
    lines = content.split("\n")

    for i, line in enumerate(lines):
        line = line.strip()
        if not line:
            continue
        if line.startswith("$ BLOCK") or line.startswith("$BLOCK"):
            if current_block:
                blocks.append(current_block)
            current_block = {"line": i+1, "header": line, "phases": [], "segments": [], "current_segment": []}
        elif line.startswith("$F0") or line.startswith("$E"):
            if current_block:
                phase_name = (
                    line.split(maxsplit=1)[1] if len(line.split()) > 1 else "Unknown"
                )
                current_block["phases"].append(phase_name)
        elif line.startswith("BLOCKEND"):
            if current_block:
                if current_block["current_segment"]:
                    current_block["segments"].append(
                        np.array(current_block["current_segment"])
                    )
                blocks.append(current_block)
                current_block = None
        elif line.startswith("BLOCK") or line.startswith("$"):
            pass
        else:
            if current_block:
                numbers = [float(x) for x in number_pattern.findall(line)]
                if len(numbers) >= 2:
                    x_val = numbers[0]
                    y_val = numbers[1]
                    if "M" in line:
                        if current_block["current_segment"]:
                            current_block["segments"].append(
                                np.array(current_block["current_segment"])
                            )
                            current_block["current_segment"] = []
                        current_block["current_segment"].append([x_val, y_val])
                    else:
                        current_block["current_segment"].append([x_val, y_val])
    if current_block:
        if current_block["current_segment"]:
            current_block["segments"].append(np.array(current_block["current_segment"]))
        blocks.append(current_block)
    return blocks

blocks = parse_thermocalc_exp(content)
print(f"Parsed {len(blocks)} data blocks.")

# =========================================================================
# 2. CURVE RECONSTRUCTION (X_SrO = 2.0 * X_overall_Sr)
# =========================================================================

# A. Miscibility gap (Block 4 contains the complete dome)
b_binodal = blocks[3] # Block 4
seg_binodal = b_binodal["segments"][0]
x_bin = seg_binodal[:, 0] * 2.0
y_bin = seg_binodal[:, 1]

# Separate into left branch (CaO-rich) and right branch (SrO-rich)
peak_idx = np.argmax(y_bin)
T_max_bin = y_bin[peak_idx]
X_peak_bin = x_bin[peak_idx]

# Left Branch (CaO-rich)
x_left = x_bin[:peak_idx + 1]
T_left = y_bin[:peak_idx + 1]
# Right Branch (SrO-rich)
x_right = x_bin[peak_idx:]
T_right = y_bin[peak_idx:]

# B. Liquidus and Solidus (High Temperature)
# Solidus: Blocks 6 and 7
sol_pts = []
for idx in [5, 6]:
    for seg in blocks[idx]["segments"]:
        for x, y in seg:
            sol_pts.append((x * 2.0, y))
sol_pts = np.array(sol_pts)
# Sort by X
sol_pts = sol_pts[np.argsort(sol_pts[:, 0])]
x_sol = sol_pts[:, 0]
y_sol = sol_pts[:, 1]

# Liquidus: Blocks 8 and 9
liq_pts = []
for idx in [7, 8]:
    for seg in blocks[idx]["segments"]:
        for x, y in seg:
            liq_pts.append((x * 2.0, y))
liq_pts = np.array(liq_pts)
# Sort by X
liq_pts = liq_pts[np.argsort(liq_pts[:, 0])]
x_liq = liq_pts[:, 0]
y_liq = liq_pts[:, 1]

# =========================================================================
# 3. EXPERIMENTAL DATA (Risold 1997 & Jacob 2000)
# =========================================================================
# Risold data et al. (1997) - AGAPT. Composition given in x_CaO -> X_SrO = 1 - x_CaO
exp_data_AGAPT = [
    (880, 0.114), (918, 0.906), (937, 0.896), (957, 0.885), (1006, 0.203),
    (1055, 0.817), (1067, 0.805), (1081, 0.790), (1097, 0.772), (1113, 0.737),
    (1117, 0.336), (1125, 0.355), (1128, 0.727)
]
exp_t_agapt = [pt[0] for pt in exp_data_AGAPT]
exp_x_agapt = [1.0 - pt[1] for pt in exp_data_AGAPT]

# Jacob data (2000) - AGAP at 1100 K. Composition in x_CaO -> X_SrO = 1 - x_CaO
exp_data_AGAP = [
    (1100, 0.760), (1100, 0.288)
]
exp_t_agap = [pt[0] for pt in exp_data_AGAP]
exp_x_agap = [1.0 - pt[1] for pt in exp_data_AGAP]

# =========================================================================
# 4. PLOT CREATION
# =========================================================================
fig, ax = plt.subplots(figsize=(8, 6), dpi=600)

# A. Plotting Calculated Curves
# Miscibility Gap (Binodal)
ax.plot(x_bin, y_bin, color="#1f77b4", linestyle="-", linewidth=1.8, label="Calculated Binodal (Halite)")

# Liquidus and Solidus
ax.plot(x_liq, y_liq, color="#d62728", linestyle="-", linewidth=1.5, label="Calculated Liquidus")
ax.plot(x_sol, y_sol, color="#e377c2", linestyle="--", linewidth=1.2, label="Calculated Solidus")

# B. Shading of Phase Regions
# 1. Miscibility gap (Halite#1 + Halite#2)
T_grid_bin = np.linspace(800, T_max_bin, 300)
# To interpolate safely, we sort the branch data by temperature
sort_l = np.argsort(T_left)
sort_r = np.argsort(T_right)
X_l_interp = np.interp(T_grid_bin, T_left[sort_l], x_left[sort_l])
X_r_interp = np.interp(T_grid_bin, T_right[sort_r], x_right[sort_r])
ax.fill_betweenx(T_grid_bin, X_l_interp, X_r_interp, color="#1f77b4", alpha=0.08)

# 2. Coexistence region L + Halite (very narrow, ~0.15 K wide)
T_grid_high = np.linspace(2870.4, 3222.5, 500)
X_sol_interp = np.interp(T_grid_high, y_sol, x_sol)
X_liq_interp = np.interp(T_grid_high, y_liq, x_liq)
ax.fill_betweenx(T_grid_high, X_sol_interp, X_liq_interp, color="#e377c2", alpha=0.3)

# C. Overlaying Experimental Points
ax.scatter(exp_x_agapt, exp_t_agapt, color="#d62728", marker="o", facecolors="none", 
           edgecolors="#d62728", s=40, zorder=5, label="AGAPT Data (Risold et al., 1997)")
ax.scatter(exp_x_agap, exp_t_agap, color="green", marker="s", s=45, zorder=5, 
           label="AGAP Data (Jacob et al., 2000)")

# D. Annotations of Singular Points and Phase Fields
# Melting points
ax.plot(0.0, 3222.5, marker="o", color="black", markersize=4)
ax.text(0.02, 3122.5, r"T$_{\mathrm{m}}$(CaO) = 3222.5 K (2949.4 °C)", fontsize=8.5, ha="left", va="center")

ax.plot(1.0, 2870.4, marker="o", color="black", markersize=4)
ax.text(0.98, 2755.0, r"T$_{\mathrm{m}}$(SrO) = 2870.4 K (2597.3 °C)", fontsize=8.5, ha="right", va="bottom")

# Miscibility Gap Critical Point
ax.plot(X_peak_bin, T_max_bin, marker="v", color="darkblue", markersize=5)
ax.text(X_peak_bin, T_max_bin + 25, rf"T$_c$ = {T_max_bin:.1f} K" + "\n" + rf"X$_{{\mathrm{{SrO}}}}$ = {X_peak_bin:.3f}", 
        fontsize=8.5, color="darkblue", ha="center", va="bottom")

# Phase field labels
ax.text(0.70, 3140, "LIQUID", fontsize=11, fontweight="bold", ha="center", color="darkred")
ax.text(0.35, 1900, "Homogeneous Solid Solution\nHalite (Ca,Sr)O", fontsize=10, ha="center", 
        bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc", alpha=0.9))
ax.text(X_peak_bin, 970, "Miscibility Gap\nHalite#1 + Halite#2", fontsize=10.5, ha="center", color="darkblue",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="#f0f8ff", edgecolor="#b0c4de", alpha=0.9))

# Annotation for the L+Halite two-phase field
ax.annotate("L + Halite", xy=(0.5, 3160.6), xytext=(0.5, 2700), fontsize=9,
            arrowprops=dict(arrowstyle="->", color="black"))

# =========================================================================
# 5. AXES CONFIGURATION (Double axes)
# =========================================================================
ax.set_xlabel(r"Mole Fraction of SrO, $X_{\mathrm{SrO}}$")
ax.set_ylabel("Temperature, $T$ (K)")
ax.set_xlim(0, 1)

# Temperature range: 600 °C to 3000 °C
T_min_K = 600.0 + 273.15
T_max_K = 3000.0 + 273.15
ax.set_ylim(T_min_K, T_max_K)

# Generate round ticks in Celsius from 600 °C to 3000 °C in steps of 200 °C
c_ticks = np.arange(600, 3001, 200)
y_ticks = c_ticks + 273.15
ax.set_yticks(y_ticks)
ax.set_yticklabels([f"{t:.0f}" for t in y_ticks])
ax.grid(True, linestyle=":", alpha=0.6, color="#cccccc")

# Configure the top axis
ax2 = ax.twiny()
ax2.set_xlim(ax.get_xlim())
ax2.set_xticks(ax.get_xticks())
ax2.set_xticklabels([f"{1.0 - t:.1f}" for t in ax.get_xticks()])
ax2.set_xlabel(r"Mole Fraction of CaO, $X_{\mathrm{CaO}}$")

# Configure the right axis
ay2 = ax.twinx()
ay2.set_ylim(ax.get_ylim())
ay2.set_yticks(y_ticks)
ay2.set_yticklabels([f"{t - 273.15:.0f}" for t in y_ticks])
ay2.set_ylabel("Temperature, $T$ (°C)")

# Legend
ax.legend(fontsize=9, loc="center left", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#cccccc")

# Title
plt.title("CaO-SrO Pseudobinary Phase Diagram\nHalite (Ca,Sr)O Miscibility Gap (600 °C to 3000 °C)", 
          fontsize=12, fontweight="bold", pad=20)

fig.tight_layout()

# Save figure
out_filename = "diagrama_CaSrO_dos_bonito.png"
out_path = out_filename
fig.savefig(out_path, bbox_inches="tight", dpi=600)
print(f"Diagram successfully saved to: {out_path}")


# =========================================================================
# 6. CREATION OF ZOOMED PLOT (2600 °C to 3000 °C)
# =========================================================================
fig2, ax_z = plt.subplots(figsize=(8, 6), dpi=600)

# A. Plotting Calculated Curves
# Liquidus and Solidus
ax_z.plot(x_liq, y_liq, color="#d62728", linestyle="-", linewidth=1.5, label="Calculated Liquidus")
ax_z.plot(x_sol, y_sol, color="#e377c2", linestyle="--", linewidth=1.2, label="Calculated Solidus")

# B. Shading of Phase Regions
# Shading L + Halite coexistence region
ax_z.fill_betweenx(T_grid_high, X_sol_interp, X_liq_interp, color="#e377c2", alpha=0.3)

# C. Annotations of Melting Points
ax_z.plot(0.0, 3222.5, marker="o", color="black", markersize=4)
ax_z.text(0.02, 3222.5, r"T$_{\mathrm{m}}$(CaO) = 3222.5 K (2949.4 °C)", fontsize=8.5, ha="left", va="center")

# Melting of SrO is at 2870.4 K (2597.3 °C) which is just below the bottom limit of 2600 °C (2873.15 K)
ax_z.plot(1.0, 2870.4, marker="o", color="black", markersize=4)

# Phase field labels
ax_z.text(0.70, 3140, "LIQUID", fontsize=11, fontweight="bold", ha="center", color="darkred")
ax_z.text(0.70, 2940, "Homogeneous Solid Solution\nHalite (Ca,Sr)O", fontsize=10, ha="center", 
          bbox=dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#cccccc", alpha=0.9))

# D. AXES CONFIGURATION (Double axes)
ax_z.set_xlabel(r"Mole Fraction of SrO, $X_{\mathrm{SrO}}$")
ax_z.set_ylabel("Temperature, $T$ (K)")
ax_z.set_xlim(0, 1)

# Zoomed Temperature range: 2600 °C to 3000 °C
T_min_z = 2600.0 + 273.15
T_max_z = 3000.0 + 273.15
ax_z.set_ylim(T_min_z, T_max_z)

# Generate ticks in Celsius from 2600 °C to 3000 °C in steps of 50 °C
c_ticks_z = np.arange(2600, 3001, 50)
y_ticks_z = c_ticks_z + 273.15
ax_z.set_yticks(y_ticks_z)
ax_z.set_yticklabels([f"{t:.0f}" for t in y_ticks_z])
ax_z.grid(True, linestyle=":", alpha=0.6, color="#cccccc")

# Configure the top axis
ax2_z = ax_z.twiny()
ax2_z.set_xlim(ax_z.get_xlim())
ax2_z.set_xticks(ax_z.get_xticks())
ax2_z.set_xticklabels([f"{1.0 - t:.1f}" for t in ax_z.get_xticks()])
ax2_z.set_xlabel(r"Mole Fraction of CaO, $X_{\mathrm{CaO}}$")

# Configure the right axis
ay2_z = ax_z.twinx()
ay2_z.set_ylim(ax_z.get_ylim())
ay2_z.set_yticks(y_ticks_z)
ay2_z.set_yticklabels([f"{t - 273.15:.0f}" for t in y_ticks_z])
ay2_z.set_ylabel("Temperature, $T$ (°C)")

# E. Legend
ax_z.legend(fontsize=9, loc="upper right", frameon=True, framealpha=0.9, facecolor="white", edgecolor="#cccccc")

# Title
plt.title("CaO-SrO Pseudobinary Phase Diagram\nZoom of Liquidus-Solidus Coexistence Region (2600 °C to 3000 °C)", 
          fontsize=12, fontweight="bold", pad=20)

# F. Inset Axis for High-Resolution Liquidus-Solidus Separation at X_SrO = 0.5
axins = ax_z.inset_axes([0.08, 0.08, 0.38, 0.38])
axins.plot(x_liq, y_liq, color="#d62728", linestyle="-", linewidth=1.5)
axins.plot(x_sol, y_sol, color="#e377c2", linestyle="--", linewidth=1.2)
axins.fill_betweenx(T_grid_high, X_sol_interp, X_liq_interp, color="#e377c2", alpha=0.3)

# Set limits for inset around X=0.5, T=3160.5 K (2887.35 °C)
axins.set_xlim(0.48, 0.52)
axins.set_ylim(3159.5, 3161.5)

# Set labels/ticks for inset
axins.set_xticks([0.48, 0.50, 0.52])
axins.set_yticks([3160, 3161])
axins.tick_params(labelsize=8)
axins.grid(True, linestyle=":", alpha=0.5)

# Add text labels to inset
axins.text(0.50, 3161.1, "LIQUID", fontsize=8, color="darkred", ha="center", fontweight="bold")
axins.text(0.50, 3159.7, "Halite", fontsize=8, color="black", ha="center", fontweight="bold")
axins.text(0.50, 3160.5, "L + Halite\n(~0.15 K width)", fontsize=7, color="purple", ha="center", va="center")

# Draw indicate_inset_zoom
ax_z.indicate_inset_zoom(axins, edgecolor="black")

fig2.tight_layout()

# Save zoomed figure
out_filename_zoom = "diagrama_CaSrO_zoom.png"
out_path_zoom = out_filename_zoom
fig2.savefig(out_path_zoom, bbox_inches="tight", dpi=600)
print(f"Zoomed diagram successfully saved to: {out_path_zoom}")

# Copy the zoomed image to the artifact folder
try:
    shutil.copy2(out_path_zoom, os.path.join(artifact_dir, out_filename_zoom))
    print(f"Copied to the artifact directory: {os.path.join(artifact_dir, out_filename_zoom)}")
except Exception as e:
    print(f"Error copying to the artifact directory: {e}")

