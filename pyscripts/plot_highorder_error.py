#-----------------------------------------------------------------------
# Python script to plot high-order advection
# schem errors and convergence rate
# Luan Santos - october 2022
#-----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from plot_scalar_field import plot

# Grid name
#gridname = 'icos_readref_sa_andes3_scvt_h1_'
gridname = 'icos_pol_scvt_h1_'

# Map projection
#map_projection = 'sphere'
map_projection = 'mercator'

# Colormap
colormap = 'seismic'

# Times
dt = ('12800','6400','3200','1600','800','400','200')

# Grid levels
grids = (1,2,3,4,5)

# errors data
errors = np.zeros((len(grids),4))

basename_ol_2 = '_HTC_trsk10_areageo_advmethodO_advorder2_tracer_error_t1036800_'
basename_ol_3 = '_HTC_trsk10_areageo_advmethodO_advorder3_tracer_error_t1036800_'
basename_ol_4 = '_HTC_trsk10_areageo_advmethodO_advorder4_tracer_error_t1036800_'
basename_gass = '_HTC_trsk10_areageo_advmethodG_tracer_error_t1036800_'

basenames = [basename_ol_2, basename_ol_3, basename_ol_4, basename_gass]
#basenames = [basename_ol_2, basename_ol_3, basename_ol_4]
#basenames = [basename_ol_2, basename_gass]

# data directory
datadir = "../data/"

# data directory
graphdir = "../graphs/"


# imodel latlon grid size
nlat = 720
nlon = 1440

for g in range(0, len(grids)):
    # Grid level
    glevel = str(g+1)

    for k in range(0,len(basenames)):
        # File to be opened
        filename = 'moist_swm_tc2_dt'+dt[g]+basenames[k]+gridname+str(glevel)+'.dat'
        print(filename)

        # Open the file and reshape
        f = open(datadir+filename,'rb')
        data = np.fromfile(f, dtype='float32')
        data = np.reshape(data,(nlat,nlon,3))

        # Get data
        val = data[:,:,2]
        errors[g,k] = np.amax(abs(val))

        plot(filename, colormap, map_projection)

# Plot the error graph
colors = ('green','blue','red','orange')
labels = ('Ollivier-Gooch 2nd order', 'Ollivier-Gooch 3rd order', 'Ollivier-Gooch 4th order', 'Gassmann')
ref = np.amax(errors)
for k in range(0, 4):
    plt.semilogy(grids, errors[:,k], color=colors[k], marker='x', label = labels[k])

# Plot reference lines
order1 = np.zeros(len(grids))
order2 = np.zeros(len(grids))
order3 = np.zeros(len(grids))
order4 = np.zeros(len(grids))

order1[0], order2[0], order3[0], order4[0] = 10.0*ref, 10.0*ref, 10.0*ref, 10.0*ref
for i in range(1, len(grids)):
    order1[i] = order1[i-1]/2.0
    order2[i] = order2[i-1]/4.0
    order3[i] = order3[i-1]/8.0
    order4[i] = order4[i-1]/16.0

plt.semilogy(grids, order1 , ':' , color='black', label = '1st order')
plt.semilogy(grids, order2 , '--', color='black', label = '2nd order')
plt.semilogy(grids, order3 , '-.', color='black', label = '3rd order')
plt.semilogy(grids, order4 , '--', color='black', label = '4rd order')
plt.xlabel('Grid level')
plt.ylabel('Error')
plt.legend()
plt.grid(True, which="both")
plt.savefig(graphdir+'errors_adv.png', format='png')
plt.close()


# Compute and plot the convergence rate
n = len(grids)
for k in range(0, 4):
    CR_linf = np.abs(np.log(errors[1:n,k])-np.log(errors[0:n-1,k]))/np.log(2.0)
    plt.ylim(0,4)
    plt.plot(grids[1:n], CR_linf, color=colors[k], marker='x', label = labels[k])
plt.xlabel('Grid level')
plt.ylabel('Convergence rate')
plt.legend()
plt.grid(True, which="both")
plt.savefig(graphdir+'convergence_rate.png', format='png')
plt.close()
