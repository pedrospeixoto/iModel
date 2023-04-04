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
run = True# Run the simulations?

# Grid levels
glevels = (4,)

# FV Schemes
mono_values = (1,) # mononotic options
fvs = ('og2','og3', 'og4', 'sg3')
rk = 'rk3'

# Grid name
#gridname = 'icos_readref_sa_andes3_scvt_h1_'
gridname = 'icos_pol_scvt_h1_'

# Plotting parameters
#map_projection = 'sphere'
map_projection = 'mercator'

# Test case - (2, 3 or 4)
TC = 3
tc = str(TC)+' 0'

# final day
if TC==2:
    fd = '12 12'
    tf = '1036800'
    dt = ('6400','3200','1000','800','400','200')
elif TC==3:
    fd = '30 30'
    tf = '2592000'
    dt = ('6400','3200','1600','800','400','200','100','50')
elif TC==5:
    fd = '8 8'
    tf = '691200'
    dt = ('1660','800','200','100','50','25')
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
fields = ('theta', 'qr', 'qv', 'qc')
field_names = (r'$\theta$', r'$q_r$', r'$q_v$', r'$q_c$')

fields_error = ('theta','qv','h','u')
field_error_name = ( r'$\theta$', r'$q_v$', r'$h$', r'$u$' )

# plot errors for different all schemes in  different norms
norm_title  = [r'$L_{\infty}$',r'$L_1$',r'$L_2$']

# error files
error_array = np.zeros((len(glevels),len(mono_values), len(fvs), len(fields_error), 2))

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

# compile the code
subprocess.run('cd .. ; make', shell=True)

for g in range(0, len(glevels)):
    # Grid level
    glevel = glevels[g]

    # update par files
    replace_line(pardir+'mesh.par', gridname+str(glevel)+'.xyz', 17)
    replace_line(pardir+'moist_swm.par', dt[glevel-1] + ' 0 0 ', 7)

    for mono in range(0, len(mono_values)):
        # update monotonic scheme
        replace_line(pardir+'moist_swm.par', str(mono_values[mono]), 25)
        for fv in range(0,len(fvs)):
            replace_line(pardir+'moist_swm.par', fvs[fv], 21)

            # Run the program
            if (run):
                subprocess.run('cd .. ; ./imodel', shell=True)

            for fd in range(0,len(fields)):
                # File to be opened
                filename = 'moist_swm_tc'+str(TC)+'_dt'+dt[glevel-1]+'_HTC_trsk10_areageo_advmethod_'+fvs[fv]
                filename_field_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(tf)+'_'+gridname+str(glevel)+'.dat'
                filename_field_t0 = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields[fd]+'_t'+str(t0)+'_'+gridname+str(glevel)+'.dat'

                # Plot the fields
                f = open(datadir+filename_field_tf,'rb')
                data_field = np.fromfile(f, dtype='float32')
                data_field = np.reshape(data_field,(nlat,nlon,3))
                val = data_field[:,:,2]
                q_min, q_max = np.amin(val), np.amax(val)
                q_min, q_max =  str("{:.2e}".format(q_min)),  str("{:.2e}".format(q_max))
                Title = field_names[fd]+' - Min = '+str(q_min)+', Max = '+str(q_max)+' - '+fvs[fv] +', mono = '+str(mono_values[mono])

                if fields[fd]=='qr' or fields[fd]=='qc':
                    q_abs = max(abs(np.amin(val)), abs(np.amax(val)))
                    q_min, q_max = -q_abs, q_abs
                    plot(filename_field_tf, 'RdBu', map_projection, qmin = q_min, qmax=q_max, title=Title)
                    #plot(filename_field_tf, 'jet', map_projection, title=Title)
                else:
                    plot(filename_field_tf, 'jet', map_projection, title=field_names[fd])

            # Plot errors - only for TC2
            if TC==2:
                for fd in range(0,len(fields_error)):
                    # File to be opened
                    filename = 'moist_swm_tc'+str(TC)+'_dt'+dt[glevel-1]+'_HTC_trsk10_areageo_advmethod_'+fvs[fv]
                    filename_error_tf = filename+'_'+rk+'_mono'+str(mono_values[mono])+'_'+fields_error[fd]+'_error_t'+str(tf)+'_'+gridname+str(glevel)+'.dat'

                    # Open the error file
                    f = open(datadir+filename_error_tf,'rb')
                    data_error = np.fromfile(f, dtype='float32')
                    data_error = np.reshape(data_error,(nlat,nlon,3))

                    # Get data
                    error_val = data_error[:,:,2]
                    #val = data_field[:,:,2]
                    #errors[g,mono,fv] = np.amax(abs(error_val))

                    # Plot
                    title = fields_error[fd]+' - '+fvs[fv]
                    eabs = max(abs(np.amin(error_val)), abs(np.amax(error_val)))
                    emin, emax = -eabs, eabs
                    plot(filename_error_tf, 'seismic', map_projection, emin, emax, title)

                    # Get errors from file
                    if fields_error[fd] != 'tracer':
                        error_file = 'moist_swm_tc'+str(TC)+'_dt'+dt[glevel-1]+'_HTC_trsk10_areageo_advmethod_'+fvs[fv]
                        error_file = error_file+'_'+rk+'_mono'+str(mono_values[mono])+'_errors_'+fields_error[fd]+'_'+gridname+str(glevel)+'.txt'
                        errors = np.loadtxt(datadir+error_file)
                        N = len(errors[:,0])
                        error_array[g,mono,fv,fd,0:2] = errors[N-1,1], errors[N-1,2]
                    else:
                        error_array[g,mono,fv,fd,0] = emax
                        error_array[g,mono,fv,fd,1] = emax
