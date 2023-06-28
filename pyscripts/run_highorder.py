#-----------------------------------------------------------------------
# Python script to run and plot high-order advection output
# Luan Santos - october 2022
#-----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from plot_scalar_field import plot
import subprocess
import os.path

# Program to be run
program = "./imodel"
run = True # Run the simulation?

# Grid levels
glevels = (2,3,4,5,6,)

# FV Schemes
fvs = ('sg3', )
#fvs = ('og2', 'og3', 'og4','sg2', 'sg3', 'sg4')
#fvs = ('og2', 'sg3')
#fvs = ('og4',)
#fvs = ('sg3',)

# mononotic filter
mono_values = (0,) 
#mono_values = (1,) 
#mono_values = (1,0) 

# time integrator
rk = 'rk3'

# Initial condition
ic = 6

# Velocity field
vf = 5

# Grid name
gridname = 'icos_pol_scvt_h1_' # uniform grid
#gridname = 'icos_readref_sa_andes3_scvt_h1_' # refined grid

# Plotting parameters
#map_projection = 'sphere'
map_projection = 'mercator'

def replace_line(filename, content, line_number):
    import re
    if os.path.exists(filename): # The file exists
        # Open the grid file
        file  = open(filename, "r")
        lines = file.readlines()

        # new content
        lines[line_number-1] = content+'\n'

        # Close the file
        file.close()

        # Write new file
        with open(filename, 'w') as file:
            for line in lines:
                file.write(line)

    else:   # The file does not exist
        print("ERROR in edit_file_line: file"+filename+" not found in /par.")
        exit()

# errors data
errors = np.zeros((len(glevels),len(mono_values),len(fvs)))

# data directory
datadir = "../data/"

# graphs directory
graphdir = "../graphs/"

# par directory
pardir = '../par/'

# imodel latlon grid size
nlat = 720
nlon = 1440

# Define high order test in mesh.par'
replace_line(pardir+'mesh.par', 'read', 5)
replace_line(pardir+'mesh.par', 'nopt', 9)
replace_line(pardir+'mesh.par', '1', 11)
replace_line(pardir+'mesh.par', '17', 15)

# Define highorder parameters
replace_line(pardir+'highorder.par', 'V', 3)
replace_line(pardir+'highorder.par', rk, 7)
replace_line(pardir+'highorder.par', 'U', 11)
replace_line(pardir+'highorder.par', str(vf), 13)
replace_line(pardir+'highorder.par', str(ic), 15)

# compile the code
#subprocess.run('cd .. ; make F90=gfortran', shell=True)
subprocess.run('cd .. ; make F90=gfortran', shell=True)

for g in range(0, len(glevels)):
    # Grid level
    glevel = glevels[g]

    # update par files
    replace_line(pardir+'mesh.par', gridname+str(glevel)+'.xyz', 17)

    for mono in range(0,len(mono_values)):
        # update monotonic scheme
        replace_line(pardir+'highorder.par', str(mono_values[mono]), 9)
        for fv in range(0,len(fvs)):
            replace_line(pardir+'highorder.par', fvs[fv], 5)
            
            string = fvs[fv]
            order = ''.join(filter(str.isdigit, string))
       
            
            
            

            # Files to be opened
            filename_phi = 'order'+order+'_v'+str(vf)+'_in'+str(ic)+'_advmethod_'+str(fvs[fv])+\
            '_'+rk+'_mono'+str(mono_values[mono])+'_phi_t_5_'+gridname+str(glevel)+'.dat'
            filename_error = 'order'+order+'_v'+str(vf)+'_in'+str(ic)+'_advmethod_'+str(fvs[fv])+\
            '_'+rk+'_mono'+str(mono_values[mono])+'_phi_t_10_'+gridname+str(glevel)+'.dat'

            # Run the program
            if (run):
                subprocess.run('cd .. ; ./imodel', shell=True)

            # Open the file and reshape
            # scalar field
            f = open(datadir+filename_phi,'rb')
            data_phi = np.fromfile(f, dtype='float32')
            data_phi = np.reshape(data_phi,(nlat,nlon,3))
            # error field
            f = open(datadir+filename_error,'rb')
            data_error = np.fromfile(f, dtype='float32')
            data_error = np.reshape(data_error,(nlat,nlon,3))

            # Get data
            error_val = data_error[:,:,2]
            val = data_phi[:,:,2]

            # store the error
            errors[g,mono,fv] = np.amax(abs(error_val))

            # Plot the fields
            q_min, q_max = np.amin(val), np.amax(val)
            q_min, q_max =  str("{:.2e}".format(q_min)),  str("{:.2e}".format(q_max))
            title = 'Min = '+str(q_min)+', Max = '+str(q_max)+', '+fvs[fv]+', mono='+str(mono_values[mono])+' \n'
            plot(filename_phi, 'jet', map_projection, -0.2, 1.2, title)

            eabs = max(abs(np.amin(error_val)), abs(np.amax(error_val)))
            emin, emax = -eabs, eabs
            plot(filename_error, 'seismic', map_projection, emin, emax)

# Plot the error graph
colors  = ('green','blue','red','orange','purple','gray')
markers = ('x','D','o','*','+','d')
labels  = fvs

emin = np.amin(abs(errors))
emax = np.amax(abs(errors))

for mono in range(0,len(mono_values)):
    for fv in range(0,len(fvs)):
        plt.semilogy(glevels, errors[:,mono,fv], color=colors[fv], marker=markers[fv], label = labels[fv])
    plt.ylim(0.95*emin, 1.05*emax)

    # Plot reference lines
    order1 = np.zeros(len(glevels))
    order2 = np.zeros(len(glevels))
    order3 = np.zeros(len(glevels))
    order4 = np.zeros(len(glevels))

    order1[0], order2[0], order3[0], order4[0] = emax, emax, emax, emax
    for i in range(1, len(glevels)):
        order1[i] = order1[i-1]/2.0
        order2[i] = order2[i-1]/4.0
        order3[i] = order3[i-1]/8.0
        order4[i] = order4[i-1]/16.0

    #plt.semilogy(glevels, order1 , ':' , color='black', label = '1st order')
    #plt.semilogy(glevels, order2 , '--', color='black', label = '2nd order')
    #plt.semilogy(glevels, order3 , '-.', color='black', label = '3rd order')
    #plt.semilogy(glevels, order4 , '--', color='black', label = '4rd order')
    plt.xticks(glevels)
    plt.xlabel('Grid level')
    plt.ylabel('Error')
    plt.legend()
    plt.grid(True, which="both")
    plt.savefig(graphdir+'errors_adv_mono'+str(mono_values[mono])+\
    '_ic'+str(ic)+'_vf'+str(vf)+'.pdf', format='pdf')
    plt.close()


    # Compute and plot the convergence rate
    n = len(glevels)
    for fv in range(0, len(fvs)):
        CR_linf = np.abs(np.log(errors[1:n,mono,fv])-np.log(errors[0:n-1,mono,fv]))/np.log(2.0)
        plt.ylim(0,4)
        plt.plot(glevels[1:n], CR_linf, color=colors[fv], marker=markers[fv], label = labels[fv])
    plt.xlabel('Grid level')
    plt.xticks(glevels[1:n])
    plt.ylabel('Convergence rate')
    plt.legend()
    plt.grid(True, which="both")
    plt.savefig(graphdir+'convergence_rate_mono'+str(mono_values[mono])+\
    '_ic'+str(ic)+'_vf'+str(vf)+'.pdf', format='pdf')
    plt.close()
