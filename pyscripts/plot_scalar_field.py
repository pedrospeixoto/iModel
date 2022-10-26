import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs

# Grid level
glevel = 5
glevel = str(glevel)

# Grid name
#gridname = 'icos_readref_sa_andes3_scvt_h1_'
#gridname = 'icos_pol_nopt_'
gridname = 'icos_pol_scvt_h1_'

# Scalar field
fieldname = 'moist_swm_tc2_dt400_HTC_trsk10_areageo_advorder3_qc_t1036800_'
#fieldname = 'moist_swm_tc4_dt200_HTC_trsk10_areageo_advorder1_qc_t2592000_'

# Map projection
#map_projection = 'sphere'
map_projection = 'mercator'

# Colormap
#colormap = 'jet'
#colormap = 'viridis'
#colormap = 'YlGnBu'
#colormap = 'seismic'
colormap = 'gist_ncar'
import matplotlib.colors as mcolors
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
             (0.125490203499794, 0.125490203499794, 0.501960813999176),
             (0.125490203499794, 0.125490203499794, 0.125490203499794),
             (0.501960813999176, 0.501960813999176, 0.501960813999176),
             (0.8784313797950745, 0.8784313797950745, 0.8784313797950745),
             (0.9333333373069763, 0.8313725590705872, 0.7372549176216125),
             (0.8549019694328308, 0.6509804129600525, 0.47058823704719543),
             (0.6274510025978088, 0.42352941632270813, 0.23529411852359772),
             (0.4000000059604645, 0.20000000298023224, 0.0)]
cmap = mcolors.ListedColormap(cmap_data, 'precipitation')
colormap = cmap

# Directories
graphdir = '../graphs/'
datadir  = '../data/'

# Figure format
fig_format = 'png'

# iModel plot latlon grid size
nlat = 720
nlon = 1440

# Filename to be opened
filename = fieldname+gridname+glevel

def plot(filename, colormap, map_projection):
    # Open the file and reshape
    f = open(datadir+filename+'.dat', 'rb')
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

plot(filename, colormap, map_projection)
