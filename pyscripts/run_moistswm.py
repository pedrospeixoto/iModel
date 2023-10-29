#-----------------------------------------------------------------------
# Python script to plot scalar field outputs from imodel
# Luan Santos - November 2022
#-----------------------------------------------------------------------

from plot_scalar_field import plot
import matplotlib.colors as mcolors
import subprocess
import os.path
import numpy as np

#----------------------------------------------------------------------------
# colorbar for clouds and precipitation
cmap_data = [(1.0, 1.0, 1.0),
             (0.3137255012989044, 0.8156862854957581, 0.8156862854957581),
             (0.0, 1.0, 1.0),
             (0.0, 0.8784313797950745, 0.501960813999176),
             (0.0, 0.7529411911964417, 0.0),
             (0.501960813999176, 0.8784313797950745, 0.0),
             (1.0, 1.0, 0.0),
             (1.0, 0.6274510025978088, 0.0),
             (1.0, 0.0, 0.0),
             (1.0, 0.125490203499794, 0.501960813999176),
             (0.9411764740943909, 0.250980406999588, 1.0),
             (0.501960813999176, 0.125490203499794, 1.0),
             (0.250980406999588, 0.250980406999588, 1.0),
             (0.125490203499794, 0.125490203499794, 0.501960813999176)]
colormap = mcolors.ListedColormap(cmap_data, 'precipitation')
#----------------------------------------------------------------------------

#----------------------------------------------------------------------------
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
#----------------------------------------------------------------------------

# Parameters
# Program to be run
program = "./imodel"
run = False # Run the simulations?
datadir = '../data/'
# Grids
#glevels = (3,3,3,3,3,3)
#glevels = (6,6,6,6,6,6)
#glevels = (7,7,7,7,7,7)
glevels = (3,3,3)
#grid_ref  = 'icos_readref_sa_andes3_scvt_h1_'
grid_unif = 'icos_pol_scvt_h1_'
gridnames=(grid_unif, grid_unif, grid_unif, grid_unif, grid_unif, grid_unif)
#gridnames=(grid_ref, grid_ref, grid_ref, grid_ref, grid_ref, grid_ref)

# FV Schemes
mono_values = (1,) # mononotic options
#fvs = ('og2','og3', 'og4','sg2', 'sg3', 'sg4')
fvs = ('sg2', 'sg3', 'sg4')
rk = 'rk3'

# Plotting parameters
#map_projection = 'sphere'
#map_projection = 'south_pole'
map_projection = 'mercator'

# Test case - (2, 3 or 4)
TC = 3
tc = str(TC)+' 0'

# final day
if TC==2:
    fd = '12 12'
    tf = '1036800'
    dt = ('800','800','800','800','400')
elif TC==3:
    fd = '30 30'
    tf = '2592000'
    #DT='3200' # ur2
    DT='1600' # ur3
    #DT='800'  # ur4 
    #DT='400'  # ur5 
    #DT='200'  # ur6 
    #DT='100'  # ur7 
    #DT='800'  # vr3
    #DT='400'  # vr4
    #DT='200'  # vr5
    #DT='100'  # vr6
    #DT='50'  # vr7
    dt = (DT,DT,DT,DT,DT,DT) #ur3

t0 = '0'

# data directory
datadir = "../data/"

# graphs directory
graphdir = "../graphs/"

# par directory
pardir = '../par/'

# imodel latlon grid size
nlat = 720
nlon = 1440

# fields to be plotted
#fields = ('theta', 'qr', 'qv', 'qc')
fields = ('qr', 'qv', 'qc')
#field_names = ('Temperature', 'Rain', 'Vapour', 'Cloud')
field_names = ('Rain', 'Vapour', 'Cloud')

fields_error = ('theta','qv','h','u')
field_error_name = ( r'$\theta$', r'$q_v$', r'$h$', r'$u$' )

# plot errors for different all schemes in  different norms
norm_title  = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']

# error files
error_array = np.zeros((len(glevels),len(mono_values), len(fields_error), 2))
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
replace_line(pardir+'moist_swm.par', '10000000000000 diam', 33)

# compile the code
subprocess.run('cd .. ; make F90=gfortran', shell=True)

for g in range(0, len(glevels)):
    # Grid level
    glevel = glevels[g]

    # update par files
    replace_line(pardir+'mesh.par', gridnames[g]+str(glevel)+'.xyz', 17)
    replace_line(pardir+'moist_swm.par', dt[g] + ' 0 0 ', 7)

    for mono in range(0, len(mono_values)):
        # update monotonic scheme
        replace_line(pardir+'moist_swm.par', str(mono_values[mono]), 25)
        replace_line(pardir+'moist_swm.par', fvs[g], 21)

        # Run the program
        if (run):
            subprocess.run('cd .. ;  export OMP_NUM_THREADS=8; ./imodel', shell=True)

        for fd in range(0,len(fields)):
            # File to be opened
            filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+dt[g]+'_HTC_trsk10_areageo_diam_hyperdiffusion_10to13.000_advmethod_'+fvs[g]
            filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(tf)+'_'+gridnames[g]+str(glevels[g])+'.dat'
            filename_field_t0 = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(t0)+'_'+gridnames[g]+str(glevels[g])+'.dat'

            # Get min/max of the fields
            f = open(datadir+filename_field_tf,'rb')
            data_field = np.fromfile(f, dtype='float32')
            data_field = np.reshape(data_field,(nlat,nlon,3))
            val = data_field[:,:,2]
            fields_min[g,mono,fd] = np.amin(val)
            fields_max[g,mono,fd] = np.amax(val)


# Plot the scalar field
for g in range(0, len(glevels)):
    for mono in range(0, len(mono_values)):
        for fd in range(0,len(fields)):
            # File to be opened
            filename = datadir+'moist_swm_tc'+str(TC)+'_dt'+dt[g]+'_HTC_trsk10_areageo_diam_hyperdiffusion_10to13.000_advmethod_'+fvs[g]
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
                plot(filename_field_tf, colormap, map_projection, qmin=Q_min, qmax=Q_max, title=Title)
                #plot(filename_field_tf, 'jet', map_projection, title=Title)
            #else:
            #    plot(filename_field_tf, 'jet', map_projection, qmin=Q_min, qmax=Q_max,  title=field_names[fd])
