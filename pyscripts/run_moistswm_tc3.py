#-----------------------------------------------------------------------
# Python script to run test case 3 of the moist shallow water model
# This is the flow over a mountain test from Zerroukat and Allen 2015
# Luan Santos - october 2023
#-----------------------------------------------------------------------

from plot_scalar_field import plot
import subprocess
import os.path
import numpy as np
from miscellaneous import *

# Parameters
# Program to be run
program = "./imodel"
run = True # Run the simulations?
#run = False # Run the simulations?

#-------------------------------------------------------------------------------------
# Grids
glevels = (3,3,3)
grid_ref  = 'icos_readref_sa_andes3_scvt_h1_'
grid_unif = 'icos_pol_scvt_h1_'
gridnames=(grid_unif, grid_unif, grid_unif, grid_ref, grid_ref, grid_ref)
gridnames=(grid_unif, grid_unif, grid_unif)

# time step
DT_unif='100' # ur7
DT_ref ='50'  # vr7
dt = (DT_unif, DT_unif, DT_unif, DT_ref, DT_ref, DT_ref)
#-------------------------------------------------------------------------------------

# FV Schemes
mono_values = (1,) # mononotic options
fvs = ('sg3', 'og3', 'og4','sg3', 'og3', 'og4')
rk = 'rk3' #time integrator

# Plotting parameters
#map_projection = 'sphere'
#map_projection = 'south_pole'
map_projection = 'mercator'

# Test case
TC = 3
tc = str(TC)+' 0'

# final day
days = 30
fd = str(days)+' '+str(days)
tf = days*86400 # seconds
t0 = '0'

# imodel latlon grid size
nlat = 720
nlon = 1440

# fields to be plotted
fields = ('qr', 'qc')
field_names = ('Rain', 'Cloud')

# min/max for plotting range
fields_min  = np.zeros((len(glevels),len(mono_values), len(fields)))
fields_max  = np.zeros((len(glevels),len(mono_values), len(fields)))

# Define high order test in mesh.par'
replace_line(pardir+'mesh.par', 'read', 5)
replace_line(pardir+'mesh.par', 'nopt', 9)
replace_line(pardir+'mesh.par', '1', 11)
replace_line(pardir+'mesh.par', '18', 15)

# Define moist swm par
replace_line(pardir+'moist_swm.par', tc, 3)
replace_line(pardir+'moist_swm.par', fd,  5)
replace_line(pardir+'moist_swm.par', 'trsk10', 11)
replace_line(pardir+'moist_swm.par', rk, 23)
replace_line(pardir+'moist_swm.par', 'geo', 27)

# Hyperdifusion parameter
replace_line(pardir+'moist_swm.par', '10000000000000 diam', 33)
hypdifname = "diam_hyperdiffusion_10to13.000"


# compile the code
subprocess.run('cd .. ; make F90=gfortran', shell=True)

for g in range(0, len(glevels)):
    # Grid level
    glevel = glevels[g]

    # update par files
    replace_line(pardir+'mesh.par', gridnames[g]+str(glevel)+'.xyz', 17)
    replace_line(pardir+'moist_swm.par', str(dt[g]) + ' 0 0 ', 7)

    for mono in range(0, len(mono_values)):
        # update monotonic scheme
        replace_line(pardir+'moist_swm.par', str(mono_values[mono]), 25)
        replace_line(pardir+'moist_swm.par', fvs[g], 21)

        # Run the program
        if (run):
            subprocess.run('cd .. ;  export OMP_NUM_THREADS=8; ./imodel', shell=True)

        for fd in range(0,len(fields)):
            # File to be opened
            filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+str(dt[g])+'_HTC_trsk10_areageo_'+hypdifname+'_advmethod_'+fvs[g]
            filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(tf)+'_'+gridnames[g]+str(glevels[g])+'.dat'
            filename_field_t0 = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(t0)+'_'+gridnames[g]+str(glevels[g])+'.dat'

            # Get min/max of the fields
            f = open(filename_field_tf,'rb')
            data_field = np.fromfile(f, dtype='float32')
            data_field = np.reshape(data_field,(nlat,nlon,3))
            val = data_field[:,:,2]
            fields_min[g,mono,fd] = np.amin(val)
            fields_max[g,mono,fd] = np.amax(val)


# Plot the scalar fields
for g in range(0, len(glevels)):
    for mono in range(0, len(mono_values)):
        for fd in range(0,len(fields)):
            # File to be opened
            filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+str(dt[g])+'_HTC_trsk10_areageo_'+hypdifname+'_advmethod_'+fvs[g]
            filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(tf)+'_'+gridnames[g]+str(glevels[g])+'.dat'
            filename_field_t0 = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(t0)+'_'+gridnames[g]+str(glevels[g])+'.dat'

            # Get min/max of the fields
            f = open(datadir+filename_field_tf,'rb')
            data_field = np.fromfile(f, dtype='float32')
            data_field = np.reshape(data_field,(nlat,nlon,3))
            val = data_field[:,:,2]
            Q_min, Q_max = np.amin(fields_min[:,:,fd]), np.amax(fields_max[:,:,fd])
            q_min, q_max = np.amin(val), np.amax(val)
            q_min, q_max =  str("{:.2e}".format(q_min)),  str("{:.2e}".format(q_max))
            Title = field_names[fd]+' - Min = '+str(q_min)+', Max = '+str(q_max)+' - '+fvs[g] +', mono = '+str(mono_values[mono])+'\n'

            if fields[fd]=='qr' or fields[fd]=='qc':
                plot(filename_field_tf, colormap2, map_projection, qmin=Q_min, qmax=Q_max, title=Title)
            else:
                plot(filename_field_tf, 'jet', map_projection, qmin=Q_min, qmax=Q_max,  title=field_names[fd])
