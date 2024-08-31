#-----------------------------------------------------------------------
# Python script to run and analyse converge of TC1 in the 
# moist shallow water model
#
# TC1 - advection of a single tracer
# IC: One gaussian hill, VF: zonal wind
# This test allows to assess the accuracy of the high order advection schemes
#
# Luan Santos - october 2023
#-----------------------------------------------------------------------

from plot_scalar_field import plot
import matplotlib.pyplot as plt
import subprocess
import os.path
import numpy as np
from miscellaneous import *

# Parameters
# Test case
TC = 1
tc = str(TC)+' 0'

# Program to be run
program = "./imodel"
#run = False # Run the simulations?
run = True # Run the simulations?

# Grids
glevels = (1,2,3,4,5,)

# Grid: unif or ref?
#grid  = 'icos_readref_sa_andes3_scvt_h1_'
grid = 'icos_pol_scvt_h1_'

# FV Schemes
mono_values = (1,) # mononotic options
#fvs = ('sg3', 'og3', 'og4')
#fvs = ('og4',)
fvs = ('sg3', 'og3', 'og4')
#fvs = ('og3', )
rk = 'rk3'

# Plotting parameters
#map_projection = 'sphere'
#map_projection = 'south_pole'
map_projection = 'mercator'


# initial time step (for grid of level 1)
if grid == 'icos_pol_scvt_h1_':
    dt0 = 6400
else:
    dt0 = 3200

# final time of integration
days = 12
fd = str(days)+' '+str(days)
tf = days*86400 # in seconds
t0 = 0

# fields to be plotted
fields = ('tracer',)
field_names = ('Tracer',)
field_errors = ('tracer_error',)
field_error_names = ('Tracer error',)


# error files
fields_min  = np.zeros((len(glevels),len(mono_values), len(fvs), len(fields)))
fields_max  = np.zeros((len(glevels),len(mono_values), len(fvs), len(fields)))
field_errors_min  = np.zeros((len(glevels),len(mono_values), len(fvs), len(field_errors)))
field_errors_max  = np.zeros((len(glevels),len(mono_values), len(fvs), len(field_errors)))
error_max  = np.zeros((len(glevels),len(mono_values), len(fvs), len(field_errors)))

# Define high order test in mesh.par'
replace_line(pardir+'mesh.par', 'read', 5)
replace_line(pardir+'mesh.par', 'nopt', 9)
replace_line(pardir+'mesh.par', '1', 11)
replace_line(pardir+'mesh.par', '18', 15)

# Define moist swm par
replace_line(pardir+'moist_swm.par', tc, 3)
replace_line(pardir+'moist_swm.par', str(fd),  5)
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

    # time step
    dt = int(dt0/2**g)

    # update par files
    replace_line(pardir+'mesh.par', grid+str(glevel)+'.xyz', 17)
    replace_line(pardir+'moist_swm.par', str(dt) + ' 0 0 ', 7)

    for fv in range(0,len(fvs)):
        for mono in range(0,len(mono_values)):
            # update monotonic scheme
            replace_line(pardir+'moist_swm.par', str(mono_values[mono]), 25)
            replace_line(pardir+'moist_swm.par', fvs[fv], 21)

            # Run the program
            if (run):
                subprocess.run('cd .. ;  export OMP_NUM_THREADS=8; ./imodel', shell=True)

            for fd in range(0,len(fields)):
                # File to be opened
                filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+str(dt)+'_HTC_trsk10_areageo_'+hypdifname+'_advmethod_'+fvs[fv]
                filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(tf)+'_'+grid+str(glevels[g])+'.dat'
                # Get min/max of the fields
                f = open(filename_field_tf,'rb')
                data_field = np.fromfile(f, dtype='float32')
                data_field = np.reshape(data_field,(nlat,nlon,3))
                val = data_field[:,:,2]
                fields_min[g,mono,fv,fd] = np.amin(val)
                fields_max[g,mono,fv,fd] = np.amax(val)

            for fd in range(0,len(field_errors)):
                # File to be opened
                filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+str(dt)+'_HTC_trsk10_areageo_'+hypdifname+'_advmethod_'+fvs[fv]
                filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+field_errors[fd]+'_t'+str(tf)+'_'+grid+str(glevels[g])+'.dat'

                # Get min/max of the fields
                f = open(filename_field_tf,'rb')
                data_field = np.fromfile(f, dtype='float32')
                data_field = np.reshape(data_field,(nlat,nlon,3))
                val = data_field[:,:,2]
                field_errors_min[g,mono,fv,fd] = np.amin(val)
                field_errors_max[g,mono,fv,fd] = np.amax(val)
                #maxabs =  max(abs(fields_min[g,mono,fv,fd]), abs(fields_max[g,mono,fv,fd]))
                #field_errors_min[g,mono,fv,fd] = -maxabs
                #field_errors_max[g,mono,fv,fd] =  maxabs
                error_max[g,mono,fv,fd] = np.amax(abs(val))

# Plot the scalar field
for g in range(0, len(glevels)):
    # time step
    dt = int(dt0/2**g)


    # scalar field plots
    for mono in range(0, len(mono_values)):
        for fv in range(0,len(fvs)):
            for fd in range(0,len(fields)):
                # File to be opened
                filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+str(dt)+'_HTC_trsk10_areageo_'+hypdifname+'_advmethod_'+fvs[fv]
                filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(tf)+'_'+grid+str(glevels[g])+'.dat'

                # Get min/max of the fields
                f = open(filename_field_tf,'rb')
                data_field = np.fromfile(f, dtype='float32')
                data_field = np.reshape(data_field,(nlat,nlon,3))
                val = data_field[:,:,2]
                Q_min, Q_max = np.amin(fields_min[:,:,:,fd]), np.amax(fields_max[:,:,:,fd])
                q_min, q_max = np.amin(val), np.amax(val)
                q_min, q_max =  str("{:.2e}".format(q_min)),  str("{:.2e}".format(q_max))
                Title = field_names[fd]+' - Min = '+str(q_min)+', Max = '+str(q_max)+' - '+fvs[fv] +', mono = '+str(mono_values[mono])+'\n'
                plot(filename_field_tf, 'jet', map_projection, qmin=Q_min, qmax=Q_max, title=Title)


    # plot errors
    for mono in range(0, len(mono_values)):
        for fv in range(0,len(fvs)):
            for fd in range(0,len(field_errors)):
                # File to be opened
                filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+str(dt)+'_HTC_trsk10_areageo_'+hypdifname+'_advmethod_'+fvs[fv]
                filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+field_errors[fd]+'_t'+str(tf)+'_'+grid+str(glevels[g])+'.dat'

                # Get min/max of the fields
                f = open(filename_field_tf,'rb')
                data_field = np.fromfile(f, dtype='float32')
                data_field = np.reshape(data_field,(nlat,nlon,3))
                val = data_field[:,:,2]
                q_min, q_max = np.amin(val), np.amax(val)
                Q_min, Q_max = np.amin(field_errors_min[g,:,:,fd]), np.amax(field_errors_max[g,:,:,fd])
                qabs = max(abs(Q_min), abs(Q_max))
                Q_min, Q_max = -qabs, qabs
                q_min, q_max =  str("{:.2e}".format(q_min)),  str("{:.2e}".format(q_max))
                Title = field_error_names[fd]+' - Min = '+str(q_min)+', Max = '+str(q_max)+' - '+fvs[fv] +', mono = '+str(mono_values[mono])+'\n'

                plot(filename_field_tf, 'seismic', map_projection, qmin=Q_min, qmax=Q_max,  title=Title)


#-------------------------------------------------------------------------------------------------
# Plot error convergence graph
colors  = ('green','blue','red','green','blue','red',)
markers = ('x','D','o','*','+','d')
linestyles = ('-','-','-','--','--','--')

for fd in range(0,len(field_errors)):
    for fv in range(0,len(fvs)):
        for mono in range(0,len(mono_values)):
            error = error_max[:,mono,fv,fd]
            plt.semilogy(glevels, error, marker='o', label = fvs[fv])
    plt.xlabel('Grid level')
    plt.ylabel('Error')
    plt.legend()
    plt.title(field_error_names[fd])
    plt.grid(True, which="both")
    plt.savefig(graphdir+'errors_'+field_errors[fd], format='pdf')
    plt.close()


# Plot convergence rate graph
for fd in range(0,len(field_errors)):
    for fv in range(0,len(fvs)):
        for mono in range(0,len(mono_values)):
            error = error_max[:,mono,fv,fd]
            n = len(error)
            CR = (np.log(error[0:n-1])-np.log(error[1:n]))/np.log(2.0)
            plt.plot(glevels[1:n], CR, marker='o', label = fvs[fv])
    plt.xlabel('Grid level')
    plt.ylabel('Convergence rate')
    plt.ylim(-0.5, 3)
    plt.legend()
    plt.title(field_error_names[fd]+' - convergence rate')
    plt.grid(True, which="both")
    plt.savefig(graphdir+'CR_'+field_errors[fd], format='pdf')
    plt.close()
