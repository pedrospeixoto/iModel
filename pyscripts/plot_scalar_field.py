#-----------------------------------------------------------------------
# Python script to plot scalar field outputs from imodel
# Luan Santos - october 2022
#-----------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

def plot(filename, colormap, map_projection):
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
    if map_projection == "mercator":
        plateCr = ccrs.PlateCarree()
        plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
    elif map_projection == "sphere":
        plateCr = ccrs.Orthographic(central_longitude=-60.0, central_latitude=0.0)
        plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi)

    plateCr._threshold = plateCr._threshold/10.
    ax = plt.axes(projection=plateCr)
    ax.stock_img()

    if map_projection == 'mercator':
        ax.gridlines(draw_labels=True)

    # Plot the scalar field
    plt.contourf(lon, lat, val, levels = 100, cmap=colormap, transform=ccrs.PlateCarree())

    # Add coastlines
    ax.coastlines()

    # Plot colorbar
    if map_projection == 'mercator':
        plt.colorbar(orientation='horizontal',fraction=0.046, pad=0.04)
    elif map_projection == 'sphere':
        plt.colorbar(orientation='vertical',fraction=0.046, pad=0.04)

    # Save the figure
    plt.savefig(graphdir+filename+'_'+map_projection+'.'+fig_format, format=fig_format)

    plt.close()
    print('Figure has been saved in '+graphdir+filename+'.'+fig_format)
