#-----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def plot(filename, colormap, map_projection, qmin=None, qmax=None, title=None):
    # Directories
    graphdir = '../graphs/'
    datadir = '../data/'

    # Figure format
    fig_format = 'png'

    # iModel plot latlon grid size
    nlat = 720
    nlon = 1440

    # Open the file and reshape
    f = open(datadir+filename, 'rb')
    data = np.fromfile(f, dtype='float32')
    data = np.reshape(data, (nlat, nlon, 3))

    # Get the data to plot
    lon = data[:,:,0]
    lat = data[:,:,1]
    val = data[:,:,2]

    # Figure quality
    dpi = 100

    # Map projection
    if map_projection == "mercator" or map_projection == 'south_america':
        plateCr = ccrs.PlateCarree()
        plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
    elif map_projection == "sphere":
        plateCr = ccrs.Orthographic(central_longitude=-60.0, central_latitude=0.0)
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi)
    elif map_projection == 'south_america':
        plateCr = ccrs.PlateCarree()
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi)

    plateCr._threshold = plateCr._threshold/10.
    ax = plt.axes(projection=plateCr)
    ax.stock_img()

    if map_projection == 'mercator' or map_projection == 'south_america':
        ax.gridlines(draw_labels=True)

    if map_projection == 'south_america':
        ax.set_extent([-100, 0, -80, -30], crs=ccrs.PlateCarree())

    if not qmin or not qmax:
        # Plot the scalar field
        plt.contourf(lon, lat, val, levels = 100, cmap=colormap, transform=ccrs.PlateCarree())
    else:
        plt.contourf(lon, lat, val, levels = np.linspace(qmin, qmax, 101), cmap=colormap, transform=ccrs.PlateCarree())

    # Plot colorbar
    if map_projection == 'mercator':
        plt.colorbar(orientation='horizontal',fraction=0.046, pad=0.04, format='%.1e')
    elif map_projection == 'sphere':
        plt.colorbar(orientation='vertical',fraction=0.046, pad=0.04, format='%.1e')
    elif map_projection == 'south_america':
        plt.colorbar(orientation='vertical',fraction=0.046, pad=0.04, format='%.1e')

    # Add coastlines
    ax.coastlines()

    # add title
    if title:
        plt.title(title)

    # Save the figure
    plt.savefig(graphdir+filename+'_'+map_projection+'.'+fig_format, format=fig_format)

    plt.close()
    print('Figure has been saved in '+graphdir+filename+'.'+fig_format)

#filename="order3_v5_in6_advmethod_sg3_rk3_mono1_phi_t_5_icos_readref_sa_andes3_scvt_h1_5.dat"
#filename="moist_swm_tc3_dt800_HTC_trsk10_areageo_advmethod_trsk_rk4_mono1_qc_t2592000_icos_pol_scvt_h1_3.dat"
#filename="moist_swm_tc3_dt800_HTC_trsk10_areageo_advmethod_trsk_rk4_mono1_qr_t2592000_icos_pol_scvt_h1_3.dat"
#filename="moist_swm_tc2_dt800_HTC_trsk10_areageo_qc_t2592000_icos_pol_scvt_h1_3.dat"
#filename="moist_swm_tc4_dt800_HTC_trsk10_areageo_qc_t2592000_icos_pol_scvt_h1_3.dat"
#filename="moist_swm_tc4_dt800_HTC_trsk10_areageo_qr_t2592000_icos_pol_scvt_h1_3.dat"
#filename = 'order2_v5_in6_mono0_phi_t_10_icos_pol_scvt_h1_4.dat'
#plot(filename, 'jet', 'mercator')
