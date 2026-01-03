import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# =========================
# Configuration parameters
# =========================

program = "./imodel"
run = False           # Run the model?
nthreads = 8          # OMP threads

glevel = 6  # Grid level
fvs = ("og2","og3","og4","sg2","sg3","sg4")              # FV schemes
mono = 0         # Monotonic filter options (0 or 1)
rk = "rk3"               # Time integrator

# Initial condition and velocity field
ic = 3
vf = 4

# Grid name
grid, gridtitle = "icos_pol_scvt_h1_", "Uniform resolution"  # uniform grid
#grid, gridtitle = "icos_readref_sa_andes3_scvt_h1_", "variable resolution (Andes)"  # refined grid

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
figformat = "pdf"

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

# =========================
# Main execution
# =========================
if __name__ == "__main__":
    import os
    import subprocess
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import cartopy.crs as ccrs
    from matplotlib.colors import BoundaryNorm

    # -------------------------
    # Compile the model
    # -------------------------
    subprocess.run("cd .. ; make F90=gfortran", shell=True)

    # -------------------------
    # FV schemes and ordering
    # -------------------------
    og_schemes = ("og2", "og3", "og4")
    sg_schemes = ("sg2", "sg3", "sg4")
    fvs_ordered = [("og2","sg2"), ("og3","sg3"), ("og4","sg4")]

    files = [
        build_filename("".join(filter(str.isdigit, fv)), vf, ic, fv, rk, mono, glevel, "phi_t_5")
        for pair in fvs_ordered for fv in pair
    ]
    titles = [fv.upper() for pair in fvs_ordered for fv in pair]

    # -------------------------
    # Run the model for each FV scheme
    # -------------------------
    for fv in [fv for pair in fvs_ordered for fv in pair]:
        update_par_files(glevel, fv, mono)
        if run:
            subprocess.run(f"cd .. ; export OMP_NUM_THREADS={nthreads}; {program}", shell=True)

    # -------------------------
    # Reference initial condition (IC)
    # -------------------------
    ic_file = build_filename("2", vf, ic, "sg2", rk, mono, glevel, "phi_t_0")
    ic_title = "Exact solution / Initial condition"

    # -------------------------
    # Create figure with GridSpec
    # -------------------------
    fig = plt.figure(figsize=(18, 22))
    gs = gridspec.GridSpec(4, 2, figure=fig, height_ratios=[1,1,1,1.2],
                           hspace=0.25, wspace=0.15)

    if grid == "icos_pol_scvt_h1_":
        lat_min, lat_max = -45, 45
        lon_min, lon_max = -80, 80
    else:
        lat_min, lat_max = -15-45, -15+45
        lon_min, lon_max = -70-80, -70+80


    # -------------------------
    # Colormap with saturation for FV methods
    # -------------------------
    cmap_sat = plt.get_cmap('jet').copy()
    cmap_sat.set_under('white')
    cmap_sat.set_over('black')
    eps = 1E-12
    levels = np.linspace(-eps, 1.01, 101)
    norm_sat = BoundaryNorm(levels, ncolors=cmap_sat.N, extend='both')

    # -------------------------
    # Plot FV schemes
    # -------------------------
    for i, fv_pair in enumerate(fvs_ordered):
        for j, fv in enumerate(fv_pair):
            ax = fig.add_subplot(gs[i,j], projection=ccrs.PlateCarree())
            with open(datadir + build_filename("".join(filter(str.isdigit, fv)), vf, ic, fv, rk, mono, glevel, "phi_t_5"), "rb") as f:
                data = np.fromfile(f, dtype='float32').reshape(nlat, nlon, 3)
            lon = data[:, :, 0]
            lat = data[:, :, 1]
            val = data[:, :, 2]
            if mono==1 :
               val[val>1.01] = 1.0
               val[val<-eps] = 0.0

            cf = ax.contourf(lon, lat, val, levels=levels, cmap=cmap_sat, norm=norm_sat,
                             extend='both', transform=ccrs.PlateCarree())
            ax.coastlines()
            vmin, vmax = val.min(), val.max()
            print(vmin,vmax)
            #ax.set_title(f"{fv.upper()}, min={vmin:.8f}, max={vmax:.3f}",
            ax.set_title(f"{fv.upper()}, min={vmin:.2e}, max={vmax:.2e}",
                         fontsize=18, pad=10)
            #ax.set_title(f"{fv.upper()}", fontsize=18, pad=10)
            ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())
            
            gl = ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
            gl.top_labels = False
            gl.right_labels = False
            gl.xlabel_style = {'size':12}
            gl.ylabel_style = {'size':12}

    # -------------------------
    # Plot reference IC spanning both columns
    # -------------------------
    ic_ax = fig.add_subplot(gs[3, :], projection=ccrs.PlateCarree())
    with open(datadir + ic_file, "rb") as f:
        data_ic = np.fromfile(f, dtype='float32').reshape(nlat, nlon, 3)
    lon_ic = data_ic[:, :, 0]
    lat_ic = data_ic[:, :, 1]
    val_ic = data_ic[:, :, 2]

    cf_ic = ic_ax.contourf(lon_ic, lat_ic, val_ic, levels=levels, cmap='jet',
                           extend='neither', transform=ccrs.PlateCarree())
    ic_ax.coastlines()
    ic_ax.set_title(ic_title, fontsize=18, pad=10)
    ic_ax.set_extent([lon_min, lon_max, lat_min, lat_max], crs=ccrs.PlateCarree())

    gl_ic = ic_ax.gridlines(draw_labels=True, linewidth=0.5, color='gray', alpha=0.7, linestyle='--')
    gl_ic.top_labels = False
    gl_ic.right_labels = False
    gl_ic.xlabel_style = {'size':14}
    gl_ic.ylabel_style = {'size':14}

    # -------------------------
    # Vertical colorbar on the right (compact)
    # -------------------------
    cbar_ax = fig.add_axes([0.935, 0.25, 0.02, 0.45])
    cbar = fig.colorbar(plt.cm.ScalarMappable(norm=norm_sat, cmap=cmap_sat),
                        cax=cbar_ax, orientation='vertical', extend='both')
    cbar.set_label(r"Tracer $\phi$ density", fontsize=14)


    # -------------------------
    # Add compact figure title (centered above panels)
    # -------------------------
    mono_str = "flux limiter ON" if mono == 1 else "flux limiter OFF"
    suptitle = (
        f"{gridtitle} (level {glevel}) | IC: {ic_name} | Wind: {vf_name} | {mono_str}"
    )

    # Fine-tune vertical position: 'y' controls the height of the title
    fig.suptitle(
        suptitle,
        fontsize=20,
        y=0.93 # << lowers the title; typical values between 0.90 and 0.94 work well
    )

    # Adjust layout to avoid overlapping with the title and the colorbar
    plt.tight_layout(rect=[0, 0, 0.93, 0.92])

    # -------------------------
    # Layout and save
    # -------------------------
    plt.tight_layout(rect=[0, 0, 0.93, 0.97])  # leave space for title and colorbar
    save_path = os.path.join(graphdir, f"{grid}level{glevel}_ic{ic}_vf{vf}_mono{mono}.png")
    plt.savefig(save_path, dpi=150, bbox_inches='tight')
    print(f"Figure saved at: {save_path}")
