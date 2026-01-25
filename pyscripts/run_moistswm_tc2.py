#-----------------------------------------------------------------------
# Python script to run and analyse converge of TC2 in the 
# moist shallow water model
#
# TC2: steady test case of Zerroukat and Allen 2015
#
# Luan Santos - october 2023
#-----------------------------------------------------------------------

import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from plot_scalar_field import plot

# =========================
# Configuration parameters
# =========================

program = "./imodel"
#run = False           # Run the model?
run = True           # Run the model?
nthreads = 14          # OMP threads

#glevels = (2, 3, 4, 5, 6, 7 )      # Grid levels
#glevels = (2, 3, 4, 5, 6 )      # Grid levels
glevels = (2, 3)      # Grid levels
#glevels = (1, 2, 3, 4, 5, 6)      # Grid levels
g0 = glevels[0]
#g0 = 7

# FV schemes
#fvs = ("sg2","sg3","sg4")
#fvs = ("og2","og3","og4")
fvs = ("og2","og3","og4","sg2","sg3","sg4")
#fvs = ("sg2","sg3",)
#fvs = ("sg2","sg2")
#fvs = ("og2","og3","og4")
#fvs = ("og2","sg2")
#fvs = ("og3","sg3")
#fvs = ("og4",)
mono_values = (2,)          # Monotonic filter options (0 or 1)

rk = "rk3"                  # Time integrator

# vector field position
HTs = ('HTC','HTC','HTC','HC','HC','HC')
HTs = ('HTC','HTC','HTC','HTC','HTC','HTC')
#HTs = ('HC','HTC')


# Grid
grid, gridtitle = "icos_pol_scvt_h1_", "uniform resolution"  # uniform grid
#grid, gridtitle = "icos_readref_sa_andes3_scvt_h1_", "variable resolution (Andes)"  # refined grid

# Define moist swm par
TC = 2
tc = str(TC)+' 0'

# final day
days = 12
fd = str(days)+' '+str(days)
tf = days*86400 # seconds
t0 = '0'



if grid == "icos_pol_scvt_h1_":
    dt0 = 5120

    # Hyperdiffusion coefficients
    hyp_type = 'diam'
    hyp_coefs_str = ('','','','','','','')
    hyp_coefs= (0,0,0,0,0,0,0,0)

else:
    dt0 = 2560

    # Hyperdiffusion coefficients
    hyp_type = 'diam'
    hyp_10to12, hyp_10to12_str = 1000000000000, hyp_type+'_hyperdiffusion_10to12.000_'
    hyp_10to13, hyp_10to13_str = 10*hyp_10to12, hyp_type+'_hyperdiffusion_10to13.000_'
    hyp_10to14, hyp_10to14_str = 10*hyp_10to13, hyp_type+'_hyperdiffusion_10to14.000_'
    hyp_10to15, hyp_10to15_str = 10*hyp_10to14, hyp_type+'_hyperdiffusion_10to15.000_'
    hyp_coefs_str = ( hyp_10to12_str, hyp_10to12_str, hyp_10to12_str, hyp_10to15_str, hyp_10to15_str, hyp_10to14_str,)
    hyp_coefs = ( hyp_10to12, hyp_10to12, hyp_10to12, hyp_10to15, hyp_10to15, hyp_10to14,)


# Test case
TC = 2
fdays = 12
tf = fdays*86400 # final time in seconds

# Fields to analyze
field_errors = ('qc', 'qv_error', 'tracer_error', 'theta_error', 'h_error')
fields = ('qc', 'qv', 'tracer', 'theta', 'h')
field_error_names = ('Cloud error', 'Vapour error', 'Tracer error', 'Temperature error', 'Fluid depth error')


# Directories
datadir = "../data/"
graphdir = "../graphs/"
pardir = "../par/"
nlat, nlon = 720, 1440
figformat = "png"
map_projection = "mercator"

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

def read_field(filename):
    """Read binary field (nlat, nlon, 3) and return the third component."""
    if not os.path.exists(filename):
        raise FileNotFoundError(f"{filename} not found")
    with open(filename, "rb") as f:
        data = np.fromfile(f, dtype="float32")
    return np.reshape(data, (nlat, nlon, 3))[:, :, 2]

def compute_norms(field):
    """Compute L_inf and L2 norms of a 2D field."""
    linf = np.max(np.abs(field))
    l2 = np.sqrt(np.mean(field**2))
    return linf, l2


def plot_error_convergence(
    glevels, errors_2d, norm_label, norm, mono, monotitle,
    gridtitle, outname, fvs, field_name, HTs
):
    colors = ("orange", "blue", "magenta", "orange", "blue", "magenta")
    markers = ("x", "D", "o", "*", "+", "d")
    linestyles = ("-", "-", "-", "--", "--", "--")

    fig, ax = plt.subplots()
    method_handles = []

    # --- Plot each FV scheme ---
    for fv_idx, fv in enumerate(fvs):
        y = errors_2d[:, fv_idx]
        HT = HTs[fv_idx]

        # Compute convergence rate (optional)
        CR = np.nan
        if len(y) > 1:
            CR = (np.log(y[-2]) - np.log(y[-1])) / np.log(2.0)

        line, = ax.semilogy(
            glevels,
            y,
            color=colors[fv_idx],
            marker=markers[fv_idx],
            linestyle=linestyles[fv_idx],
            label=f"{fv}" if not np.isnan(CR) else fv
        )
        method_handles.append(line)

    # --- Reference lines (1st and 2nd order) ---
    gsub = glevels[-2:]
    emax = 10.0 * np.max(errors_2d[-2:, :])
    orders = {"1st order": ("-.", 2), "2nd order": ("-", 4)}
    ref_handles = []
    for label, (style, factor) in orders.items():
        ref = [emax, emax / factor]
        line, = ax.semilogy(gsub, ref, style, color="black", label=label)
        ref_handles.append(line)

    # --- Axis settings ---
    ax.set_xticks(glevels)
    ax.set_xlabel("Grid level", fontsize=14)
    ax.set_ylabel(f"{norm_label} error", fontsize=14)
    ax.set_title(f"{field_name} - {monotitle} - {gridtitle}", fontsize=12)
    ax.grid(True, which="both")
    ax.tick_params(axis='x', labelsize=14)
    ax.tick_params(axis='y', labelsize=14)


    # --- Two legends ---
    # 1) FV schemes — top-right corner
    leg1 = ax.legend(
        handles=method_handles,
        fontsize=10,
        loc="upper right",
        bbox_to_anchor=(0.98, 0.98),
        frameon=True
    )
    # 2) Reference lines — bottom-left corner
    leg2 = ax.legend(
        handles=ref_handles,
        fontsize=10,
        loc="lower left",
        bbox_to_anchor=(0.02, 0.02),
        frameon=True
    )
    ax.add_artist(leg1)  # keep both legends visible

    plt.tight_layout()
    plt.savefig(graphdir + outname + f"_mono{mono}_norm{norm}.{figformat}")
    plt.close(fig)


def plot_convergence_rate(glevels, errors_2d, norm_label, norm, mono, monotitle, gridtitle, outname, fvs, field_name, HTs):
    colors = ("orange", "blue", "magenta", "orange", "blue", "magenta")
    markers = ("x", "D", "o", "*", "+", "d")
    linestyles = ("-", "-", "-", "--", "--", "--")

    fig, ax = plt.subplots()
    for fv_idx, fv in enumerate(fvs):
        HT = HTs[fv_idx]
        err = errors_2d[:, fv_idx]
        with np.errstate(divide='ignore', invalid='ignore'):
            CR = (np.log(err[:-1]) - np.log(err[1:])) / np.log(2.0)
        ax.plot(glevels[1:], CR, color=colors[fv_idx], marker=markers[fv_idx],
                linestyle=linestyles[fv_idx], label=fv)

    ax.set_xticks(glevels[1:])
    ax.set_xlabel("Grid level", fontsize=14)
    ax.set_ylabel(f"{norm_label} convergence rate", fontsize=14)
    ax.set_title(f"{field_name} - {monotitle} - {gridtitle} ", fontsize=12)
    ax.grid(True, which="both")
    #ax.set_ylim(-0.5, 3)
    ax.legend(fontsize=11)
    plt.savefig(graphdir + "cr_" + outname + f"_mono{mono}_norm{norm}.{figformat}")
    plt.close(fig)

# =========================
# Main execution
# =========================

# Compile model
subprocess.run("cd .. ; make F90=gfortran", shell=True)

shape = (len(glevels), len(fvs), len(mono_values), len(field_errors))
errors_linf = np.zeros(shape)
errors_l2 = np.zeros(shape)

# Set moist swm case in mesh.par
replace_line(pardir + "mesh.par", "18", 15)

# Define some moist swm parameters
replace_line(pardir+'moist_swm.par', tc, 3)
replace_line(pardir+'moist_swm.par', fd,  5)

# Shallow-water solver parameters
cell_vec_recon = 'trsk'
coriolis = 'trsk'
scalar_interp = 'trsk'
grad = 'trsk'
replace_line(pardir+'moist_swm.par', 'none', 11)            # Wrapper
replace_line(pardir+'moist_swm.par', f'{cell_vec_recon} 1.0', 13)        # Cell vec reconstruction method (Kenergy) / Gassmann parameter
replace_line(pardir+'moist_swm.par', coriolis, 15)          # Coriolis vec reconstruction method
replace_line(pardir+'moist_swm.par', scalar_interp, 17)     # Scalar interpolations
replace_line(pardir+'moist_swm.par', grad, 19)              # Gradient discrete method

for g_idx, glevel in enumerate(glevels):
    dt = dt0 // 2**(glevel-1)

    # hyperdiffusion coefficients
    hyper_diff_coef = hyp_coefs[g_idx]
    hyper_diff_coef_str = hyp_coefs_str[g_idx]

    # update mesh and parameters
    replace_line(pardir+'mesh.par', grid+str(glevel)+'.xyz', 17)
    replace_line(pardir+'moist_swm.par', f"{dt} 0 0", 7)
    replace_line(pardir+'moist_swm.par', f"{hyper_diff_coef} {hyp_type}", 33)

    for fv_idx, fv in enumerate(fvs):
        HT = HTs[fv_idx]
        replace_line(pardir+'moist_swm.par', f"{HT}", 9)
        for m_idx, mono in enumerate(mono_values):
            replace_line(pardir+'moist_swm.par', str(mono), 25)
            replace_line(pardir+'moist_swm.par', fv, 21)

            # Run model
            if run and glevel >= g0:
                subprocess.run(f"cd .. ; export OMP_NUM_THREADS={nthreads}; {program}", shell=True)

            for fd_idx, (fd_error, fd) in enumerate(zip(field_errors, fields)):
                swm_pars = f"vrec{coriolis}_crec{cell_vec_recon}_sint{scalar_interp}_grd{grad}"
                filename = f"{datadir}moist_swm_tc{TC}_dt{dt}_{HT}_{swm_pars}_areageo_{hyper_diff_coef_str}advmethod_{fv}_{rk}_mono{mono}_{fd}_t{tf}_{grid}{glevel}.dat"
                filename_error = f"{datadir}moist_swm_tc{TC}_dt{dt}_{HT}_{swm_pars}_areageo_{hyper_diff_coef_str}advmethod_{fv}_{rk}_mono{mono}_{fd_error}_t{tf}_{grid}{glevel}.dat"

                error_val = read_field(filename_error)

                error_linf, error_l2 = compute_norms(error_val)
                if (fd=='qc'): linf, l2 = 1.0, 1.0
                else:
                    val = read_field(filename)
                    linf, l2 = compute_norms(val)

                errors_linf[g_idx, fv_idx, m_idx, fd_idx] = error_linf/linf
                errors_l2[g_idx, fv_idx, m_idx, fd_idx] = error_l2/l2

# ---------------------------
# Plot convergence for each field and mono
for fd_idx, fd in enumerate(field_errors):
    for m_idx, mono in enumerate(mono_values):
        errors_linf_slice = errors_linf[:, :, m_idx, fd_idx]
        errors_l2_slice = errors_l2[:, :, m_idx, fd_idx]
        outname = f"{fd}_{grid}_{fd}_TC{TC}_{rk}"
        if mono==0:
            monotitle='Without flux limiter'
        else:
            monotitle='With flux limiter'

        title = f"{field_error_names[fd_idx]} - {monotitle} - {gridtitle}"
        plot_error_convergence(glevels, errors_linf_slice, r"$L_\infty$", 'linf', mono, monotitle, gridtitle, outname, fvs, field_error_names[fd_idx], HTs)
        plot_error_convergence(glevels, errors_l2_slice, r"$L_2$", 'l2', mono, monotitle, gridtitle, outname, fvs, field_error_names[fd_idx], HTs)
        plot_convergence_rate(glevels, errors_linf_slice, r"$L_\infty$", 'linf', mono, monotitle, gridtitle, outname, fvs, field_error_names[fd_idx], HTs)
        plot_convergence_rate(glevels, errors_l2_slice, r"$L_2$", 'l2', mono, monotitle, gridtitle, outname, fvs, field_error_names[fd_idx], HTs)

print("[INFO] Done.")

