import os
import numpy as np
import matplotlib.pyplot as plt

# -----------------------------
# Configuration
# -----------------------------
advmethods = ["og2", "og3", "og4", "sg2", "sg3", "sg4"]  # advection methods
glevels = [2, 3, 4, 5, 6, 7]                                      # grid levels to read
datadir = "../data/"                                     # data directory
graphdir = "../graphs/"                                  # output directory

# create graph directory if it does not exist
os.makedirs(graphdir, exist_ok=True)

# Initial condition and velocity field
ic = 6
vf = 5
mono = 1  # 1 = with flux limiter, 0 = without
grid, gridtitle = "icos_pol_scvt_h1_", "uniform resolution"  # uniform grid
grid, gridtitle = "icos_readref_sa_andes3_scvt_h1_", "variable resolution (Andes)"  # refined grid
# -----------------------------
# Read data
# -----------------------------
# Structure: data[glevel][advmethod] = Nx2 array (time, error)
data = {}

for g in glevels:
    data[g] = {}
    for method in advmethods:
        order = f"order{method[-1]}"  # extract order from method name
        filename = f"{order}_v{vf}_in{ic}_advmethod_{method}_rk3_mono{mono}_{grid}{g}_errors_evol.txt"
        filepath = os.path.join(datadir, filename)
        if os.path.isfile(filepath):
            data[g][method] = np.loadtxt(filepath)
        else:
            print(f"Warning: file {filepath} not found. Skipping.")

# -----------------------------
# Plotting style configuration
# -----------------------------
colors = ("orange", "blue", "magenta", "orange", "blue", "magenta")
linestyles = ("-", "-", "-", "--", "--", "--")

# -----------------------------
# Plot temporal error for all grid levels and save
# -----------------------------
mono_label = "with flux limiter" if mono == 1 else "without flux limiter"

for glevel in glevels:
    plt.figure(figsize=(8,6))
    for i, method in enumerate(advmethods):
        if method in data[glevel]:
            step = 10*glevel
            t = data[glevel][method][:,0]
            err = data[glevel][method][:,1]
            plt.plot(t[::step], err[::step], color=colors[i], linestyle=linestyles[i], label=method)

    plt.ylim(10**(-5), 10**1)
    plt.xlabel("Time", fontsize=14)
    plt.ylabel("Error", fontsize=14)
    plt.title(f"Temporal error evolution for {gridtitle} with glevel={glevel}\nGaussian hill, zonal wind, {mono_label}", fontsize=14)
    plt.yscale("log")
    plt.grid(True, which="both", linestyle="--", linewidth=0.5)
    plt.legend(fontsize=20)
    plt.tight_layout()

    # Increase tick font size
    plt.tick_params(axis='both', which='major', labelsize=16)  # major ticks
    plt.tick_params(axis='both', which='minor', labelsize=16)  # minor ticks (if any)

    # save figure with IC, VF, mono, grid and glevel in filename
    save_filename = f"error_ic{ic}_vf{vf}_mono{mono}_{grid}_glevel{glevel}.png"
    savepath = os.path.join(graphdir, save_filename)
    plt.savefig(savepath, dpi=300)
    plt.close()

