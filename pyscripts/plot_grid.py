####################################################################################
#
# Plotting routines.
#
# Luan da Fonseca Santos - July 2022
# (luan.santos@usp.br)
####################################################################################

import cartopy.crs as ccrs
import numpy as np
import matplotlib.pyplot as plt

####################################################################################
# This routine plots the icos grid.
####################################################################################
fig_format = 'png' # Figure format
def plot_grid(map_projection, lats, lons, N, name):
   # Figure resolution
   dpi = 100

   # Grid directory
   graphsdir = "../graphs/"

   print("--------------------------------------------------------")
   print('Plotting latlon grid using '+map_projection+' projection...')
   if map_projection == "mercator":
      plateCr = ccrs.PlateCarree()
      plt.figure(figsize=(1832/dpi, 977/dpi), dpi=dpi)
   elif map_projection == 'sphere':
      plateCr = ccrs.Orthographic(central_longitude=-60.0, central_latitude=40.0)
      plt.figure(figsize=(800/dpi, 800/dpi), dpi=dpi)
   else:
      print('ERROR: Invalid map projection.')
      exit()

   plateCr._threshold = plateCr._threshold/10.
   ax = plt.axes(projection=plateCr)
   ax.stock_img()

   for i in range(0, N, 2):
       # Plot vertices
       A_lon, A_lat = lons[i],   lats[i]
       B_lon, B_lat = lons[i+1], lats[i+1]

       plt.plot([A_lon, B_lon], [A_lat, B_lat],linewidth=1, color='blue', transform=ccrs.Geodetic())

   #if map_projection == 'mercator':
   #   ax.gridlines(draw_labels=True)

   # Get quadrature nodes
   quadrature_nodes = np.loadtxt("../grid/icos_pol_scvt_h1_2_gaussnodes.gmt")
   quadlons = quadrature_nodes[:,0]
   quadlats = quadrature_nodes[:,1]
   Nquad = len(quadlons)

   for i in range(0, Nquad):
       # Plot vertices
       A_lon, A_lat = quadlons[i], quadlats[i]
       #print(A_lon, A_lat)
       plt.plot(A_lon, A_lat, marker='+', color='red')

   #for i in range(120, Nquad):
       # Plot vertices
   #    A_lon, A_lat = quadlons[i], quadlats[i]
       #print(A_lon, A_lat)
   #    plt.plot(A_lon, A_lat, marker='X', color='gray')



   # Save the figure
   plt.savefig(graphsdir+name+"_"+map_projection+'.'+fig_format, format=fig_format)
   print('Figure has been saved in latlon_'+map_projection+'.'+fig_format)
   print("--------------------------------------------------------\n")
   plt.close()


# Grid directory
griddir = "../grid/"

# Map projection
#map_projection = "sphere"
map_projection = "mercator"

# Grid name (containing edges info)
name = "icos_pol_scvt_h1_2_edhx"
num_lines = sum(1 for line in open(griddir+name+".gmt"))

# Open grid file
f = open(griddir+name+".gmt",'r')
f.readline()

# Number of edges
n_edges = (num_lines-1)//3
grid_vertices = np.zeros((2*n_edges,2))

k = 0
for i in range(0, n_edges):
    # Read a line
    f.readline()

    # Edge starting point
    line = f.readline()

    # Get vertice longitude
    grid_vertices[k,0] = float(line[0:18])

    # Get vertices latitude
    grid_vertices[k,1] = float(line[18:])

    # Edge  endpoint
    k = k+1
    line = f.readline()

    # Get vertice longitude
    grid_vertices[k,0] = float(line[0:18])

    # Get vertices latitude
    grid_vertices[k,1] = float(line[18:])

    k = k+1
   #print(f.readline())

lons = grid_vertices[:,0]
lats = grid_vertices[:,1]
N = len(lats)

plot_grid(map_projection, lats, lons, N, name)
