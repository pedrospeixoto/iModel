# -----------------------------------------------------------------------
# Python script to run and analyse convergence of high-order advection
# schemes in the imodel framework
#
# Example: advection of two Gaussian hills under a divergent flow
#
# Luan Santos - Rewritten in modern style, October 2025
# Adjusted to include mono dimension in shape
# -----------------------------------------------------------------------

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from plot_scalar_field import plot

# =========================
# Configuration parameters
# =========================

program = "./imodel"
run = True           # Run the model?
nthreads = 14          # OMP threads

glevels = (2, 3, 4, 5, 6, 7)      # Grid levels
fvs = ("og2","og3","og4","sg2","sg3","sg4")              # FV schemes
mono_values = (0, 1,)          # Monotonic filter options (0 or 1)
rk = "rk3"                  # Time integrator
g0 = glevels[0]
#g0 = 4


# Initial condition and velocity field
ic = 2
vf = 4

# Grid name
grid, gridtitle = "icos_pol_scvt_h1_", "uniform resolution"  # uniform grid
grid, gridtitle = "icos_readref_sa_andes3_scvt_h1_", "variable resolution (Andes)"  # refined grid

# Derived names
if ic == 2:
    ic_name = "Two Gaussian hills"
elif ic == 3:
    ic_name = "Slotted cylinder"
elif ic == 6:
    ic_name = "One Gaussian hill"
else:
    ic_name = ""

if vf == 3:
    vf_name = "divergent deformational flow"
elif vf == 4:
    vf_name = "div.-free deform. + zonal wind"
elif vf == 5:
    vf_name = "zonal wind"
else:
    vf_name = ""

testname = f"{ic_name} - {vf_name}"

# =========================
# Directories and constants
# =========================
datadir = "../data/"
graphdir = "../graphs/"
pardir = "../par/"
nlat, nlon = 720, 1440
map_projection = "mercator"
figformat = "png"

# =========================
# Helper functions
# =========================

def replace_line(filename, content, line_number):
    """Replace a specific line in a file with new content."""
    if not os.path.exists(filename):
        print(f"ERROR: file {filename} not found.")
        exit(1)
    with open(filename, "r") as f:
        lines = f.readlines()
    lines[line_number - 1] = content + "\n"
    with open(filename, "w") as f:
        f.writelines(lines)


def update_par_files(glevel, fv, mono):
    """Update parameter files for a given configuration."""
    # mesh.par
    replace_line(pardir + "mesh.par", "read", 5)
    replace_line(pardir + "mesh.par", "nopt", 9)
    replace_line(pardir + "mesh.par", "1", 11)
    replace_line(pardir + "mesh.par", "17", 15)
    replace_line(pardir + "mesh.par", f"{grid}{glevel}.xyz", 17)

    # highorder.par
    replace_line(pardir + "highorder.par", "V", 3)
    replace_line(pardir + "highorder.par", fv, 5)
    replace_line(pardir + "highorder.par", str(mono), 9)
    replace_line(pardir + "highorder.par", rk, 7)
    replace_line(pardir + "highorder.par", "U", 11)
    replace_line(pardir + "highorder.par", str(vf), 13)
    replace_line(pardir + "highorder.par", str(ic), 15)


def build_filename(order, vf, ic, fv, rk, mono, glevel, suffix):
    """Construct standard filename for data files."""
    return (
        f"order{order}_v{vf}_in{ic}_advmethod_{fv}_{rk}_mono{mono}_{suffix}_"
        f"{grid}{glevel}.dat"
    )


def read_field(filename):
    """Read binary field (nlat, nlon, 3) and return the third component."""
    with open(datadir + filename, "rb") as f:
        data = np.fromfile(f, dtype="float32")
    return np.reshape(data, (nlat, nlon, 3))[:, :, 2]


def compute_norms_from_field(field):
    """Compute L_inf and L2 norms (simple domain-averaged L2)."""
    linf = np.max(np.abs(field))
    l2 = np.sqrt(np.mean(field**2))
    return linf, l2


def plot_error_convergence(glevels, errors_2d, norm_label, norm, mono, monotitle, gridtitle, outname):
    """
    Plot error vs grid level with reference lines and two separate legends.

    Parameters
    ----------
    glevels : tuple/list
        Grid levels.
    errors_2d : ndarray
        Error values with shape (len(glevels), len(fvs)).
    norm_label : str
        Label for the norm, e.g., "L_inf" or "L2".
    norm : str
        Norm identifier (used in filename).
    mono : int
        Monotonic filter value.
    outname : str
        Base name for output figure.
    """
    import matplotlib.pyplot as plt
    import numpy as np

    # ----------------------------
    # Plotting style configuration
    # ----------------------------
    colors = ("orange", "blue", "magenta", "orange", "blue", "magenta")
    markers = ("x", "D", "o", "*", "+", "d")
    linestyles = ("-", "-", "-", "--", "--", "--")

    if len(glevels) < 2:
        raise ValueError("Need at least two grid levels to plot convergence.")

    fig, ax = plt.subplots()

    method_handles = []  # handles for the first legend (methods)
    ref_handles = []     # handles for the reference lines

    # ----------------------------
    # Plot errors for each FV scheme
    # ----------------------------
    for fv_idx, fv in enumerate(fvs):
        y = errors_2d[:, fv_idx]
        mask = y > 0
        CR = np.nan
        if len(mask) > 1:
            try:
                # Compute convergence rate using last two points
                CR = (np.log(errors_2d[-2, fv_idx]) - np.log(errors_2d[-1, fv_idx])) / np.log(2.0)
            except Exception:
                pass

        line, = ax.semilogy(
            glevels,
            errors_2d[:, fv_idx],
            color=colors[fv_idx],
            marker=markers[fv_idx],
            linestyle=linestyles[fv_idx],
            label=f"{fvs[fv_idx]} - order {CR:.1f}" if not np.isnan(CR) else f"{fvs[fv_idx]}"
        )
        method_handles.append(line)

    # ----------------------------
    # Add reference lines (2nd to 4th order)
    # ----------------------------
    gsub = glevels[-2:]               # take last two grid levels
    emax = np.max(errors_2d[-2:, :])  # scale based on largest error
    emin = np.min(errors_2d[-2:, :])  # scale based on largest error
    orders = {
        "1st order": ("--", 2),
        "2nd order": ("-.", 4),
        "3rd order": ("-", 8),
    }
    #"4th order": (":", 16),

    for label, (style, factor) in orders.items():
        ref = [emax, emax / factor]
        #ref = [emin, emin / factor]
        line, = ax.semilogy(gsub, ref, style, color="black", label=label)
        ref_handles.append(line)

    # ----------------------------
    # Configure axes
    # ----------------------------
    ax.set_xticks(glevels)
    ax.xaxis.set_major_locator(plt.MaxNLocator(integer=True))
    ax.set_xlabel("Grid level", fontsize=14)
    ax.set_ylabel(f"{norm_label} error", fontsize=14)
    ax.set_title(f"{testname}\n{monotitle} - {gridtitle}", fontsize=12)
    ax.grid(True, which="both")
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)

    # ----------------------------
    # Add legends
    # ----------------------------
    # First legend: FV schemes and estimated order
    leg1 = ax.legend(handles=method_handles, fontsize=11, loc='lower left')#, title="Methods")
    ax.add_artist(leg1)

    # Second legend: reference lines
    ax.legend(handles=ref_handles, fontsize=11, loc='lower center')

    # ----------------------------
    # Save figure
    # ----------------------------
    plt.savefig(graphdir + outname + f"_mono{mono}_norm{norm}." + figformat)
    plt.close(fig)



def plot_convergence_rate(glevels, errors_2d, norm_label, norm, mono, monotitle, gridtitle, outname):
    """
    Plot convergence rate.
    errors_2d has shape (len(glevels), len(fvs)).
    """
    colors = ("orange", "blue", "magenta", "orange", "blue", "magenta")
    markers = ("x", "D", "o", "*", "+", "d")
    linestyles = ("-", "-", "-", "--", "--", "--")

    for fv_idx, fv in enumerate(fvs):
        # compute convergence rates between consecutive grid levels
        err = errors_2d[:, fv_idx]
        # avoid zeros or negative/inf values
        with np.errstate(divide='ignore', invalid='ignore'):
            CR = (np.log(err[:-1]) - np.log(err[1:])) / np.log(2.0)
        plt.plot(
            glevels[1:], CR, color=colors[fv_idx], marker=markers[fv_idx],
            linestyle=linestyles[fv_idx], label=fvs[fv_idx]
        )

    plt.xticks(glevels[1:])
    plt.xlabel("Grid level", fontsize=14)
    plt.ylabel(f"{norm_label} convergence rate", fontsize=14)
    plt.title(f"{testname}\n{monotitle} - {gridtitle}", fontsize=12)
    plt.grid(True, which="both")
    plt.legend(fontsize=11)
    plt.savefig(graphdir + "cr_" + outname + f"_mono{mono}_norm{norm}." + figformat)
    plt.close()


# =========================
# Main execution
# =========================
if __name__ == "__main__":

    # Compile
    subprocess.run("cd .. ; make F90=gfortran", shell=True)

    # shape now includes glevels, fvs (schemes) and monos
    shape = (len(glevels), len(fvs), len(mono_values))
    print(f"[INFO] error arrays shape (glevels, fvs, monos): {shape}")

    # allocate error arrays with 3 dimensions
    errors_linf = np.zeros(shape)
    errors_l2 = np.zeros(shape)

    for g_idx, glevel in enumerate(glevels):
        for fv_idx, fv in enumerate(fvs):
            for m_idx, mono in enumerate(mono_values):

                update_par_files(glevel, fv, mono)
                order = "".join(filter(str.isdigit, fv))

                # Run model
                if run and glevel>=g0:
                    subprocess.run(
                        f"cd .. ; export OMP_NUM_THREADS={nthreads}; {program}",
                        shell=True,
                    )

                # Filenames
                fname_phi = build_filename(order, vf, ic, fv, rk, mono, glevel, "phi_t_5")
                fname_error = build_filename(order, vf, ic, fv, rk, mono, glevel, "phi_t_10")
                fname_error_norms = (
                    f"order{order}_v{vf}_in{ic}_advmethod_{fv}_{rk}_mono{mono}_"
                    f"{grid}{glevel}_errors.txt"
                )

                # Read data (these will raise if files are missing — keep original behavior)
                val = read_field(fname_phi)
                err = read_field(fname_error)
                enorms = np.loadtxt(datadir + fname_error_norms)

                # Store norms into 3D arrays
                #print(datadir + fname_error_norms)
                #print(enorms)
                #exit()
                errors_linf[g_idx, fv_idx, m_idx] = enorms[0]
                errors_l2[g_idx, fv_idx, m_idx] = enorms[1]

                # Plot scalar field and error for each case (same behavior as original)
                #plot(datadir + fname_phi, "jet", map_projection, -0.2, 1.2,
                #     f"FV={fv}, mono={mono}\nMin={val.min():.2e}, Max={val.max():.2e}")
                eabs = np.max(np.abs(err))
                #plot(datadir + fname_error, "seismic", map_projection, -eabs, eabs,
                #     f"Error - FV={fv}, mono={mono}")

    # ---------------------------
    # Convergence plots (now per mono)
    # ---------------------------
    for m_idx, mono in enumerate(mono_values):
        # slice arrays for this mono: shape -> (len(glevels), len(fvs))
        errors_linf_slice = errors_linf[:, :, m_idx]
        errors_l2_slice = errors_l2[:, :, m_idx]

        outname = f"adv_{grid}ic{ic}_vf{vf}_{rk}"

        # print the plotting grid shape (glevels x fvs) for user's info
        plot_shape_2d = (len(glevels), len(fvs))
        print(f"[INFO] plotting convergence for mono={mono}, 2D shape: {plot_shape_2d}")
        if mono == 1:
            monotitle = 'With flux limiter'
        else:
            monotitle = 'Without flux limiter'

        plot_error_convergence(glevels, errors_linf_slice, r"$L_\infty$", 'linf', mono, monotitle, gridtitle, outname)
        plot_error_convergence(glevels, errors_l2_slice, r"$L_2$", 'l2', mono, monotitle, gridtitle, outname)
        plot_convergence_rate(glevels, errors_linf_slice, r"$L_\infty$", 'linf', mono, monotitle, gridtitle, outname)
        plot_convergence_rate(glevels, errors_l2_slice, r"$L_2$", 'l2', mono, monotitle, gridtitle, outname)

    print("[INFO] Done.")

