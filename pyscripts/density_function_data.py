#! /usr/bin/env python3
#----------------------------------------------------
# This script generates the density function data
# defined in a lat-lon grid
# Output: densf_table.dat (in altitude directory)
# To use this data in iModel grids generator, you must
# set the following options in the mesh.par file:
#
# !Kind of mesh use: icos, octg, read, rand
# icos
# !Position: eqs, pol, ran, ref, readref, readref_andes
# readref
# !Optimization of mesh: nopt, sprg, scvt, salt, hr95
# scvt
#----------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as ss

#----------------------------------
#converts degrees to radians
deg2rad = np.pi/180.0
#----------------------------------
#converts degrees to radians
rad2deg = 1.0/deg2rad
#----------------------------------

#-----------------------------------------------------------------------
#Transforms geographical coordinates (lat,lon) to Cartesian coordinates.
#-----------------------------------------------------------------------
def sph2cart(lat,lon):
    coslat = np.cos(lat)
    x = coslat*np.cos(lon)
    y = coslat*np.sin(lon)
    z = np.sin(lat)
    return x,y,z

#-----------------------------------------------------------------------
# Compute the density function
#-----------------------------------------------------------------------
def density_function(lat,lon):
   #Parameters
   #Center of circular refinement region in lat-lon
   latc = 0.0*deg2rad
   lonc = 0.0*deg2rad
      
   #Center of circular refinement region in R^3 coordinates
   coslat = np.cos(latc)
   xc = coslat*np.cos(lonc)
   yc = coslat*np.sin(lonc)
   zc = np.sin(latc)

   #radius of circular refinement region
   maxdist = 10.0*deg2rad

   #width of transition zone
   epsilons = 25.0*deg2rad
   
   #gammas is approximately the ratio between high and low resolution cell diameters   
   gammas = 10.0
   
   #convert lat-lon points to (x,y,z) coordinates
   [x,y,z] = sph2cart(lat,lon)
   
   #distance to center metric (Euclidian norm in R^3)
   dists = np.sqrt((x-xc)**2 + (y-yc)**2 + (z-zc)**2)
   
   #auxilary function
   sx = (dists-maxdist)/epsilons

   #Boolean indexes
   #point closer to the center
   closer_to_center = (dists <= maxdist) 

   #point in transition zone
   transition_zone = (dists <= maxdist + epsilons) & (dists > maxdist)
   
   #point far from the center
   far_from_center = (dists > maxdist + epsilons)
   
   #set density
   dens_f = np.zeros(np.shape(lat))
   dens_f[closer_to_center] = gammas**4
   dens_f[transition_zone]  = ((1.0-sx[transition_zone])*gammas+sx[transition_zone])**4
   dens_f[far_from_center]  = 1.0
 
   #normalization - make it in [0,1]
   dens_f = dens_f/gammas**4
   return dens_f

#lat-lon grid
nlat = 1024
nlon = 2048

lats = np.linspace(-90.0,90.0,nlat+1)
lons = np.linspace(-180.0,180.0,nlon+1)

#converts degrees to radians
deg2rad = np.pi/180.0
lats = lats*deg2rad
lons = lons*deg2rad
lon, lat = np.meshgrid(lons, lats)

#compute the density function in the lat-lon grid
print("Computing the density function in the lat-lon grid...")
data = density_function(lat,lon)
print("Done\n")

#saves the data in densf_table
datadir = "../altitude/"
print("Saving the data in the file "+datadir+"densf_table.dat...")
with open(datadir+'densf_table.dat', 'w') as f:
    f.write(str(nlat+1))  
    f.write(" ")                
    f.write(str(nlon+1))                                
    f.write('\n')  
    for i in range(0,nlat+1):
        for j in range(0,nlon+1):
                f.write(str(lat[i,j]))
                f.write(" ")                
                f.write(str(lon[i,j]))  
                f.write(" ")                
                f.write(str(data[i,j]))                                
                f.write('\n')                
print("Done\n")

#plot the data
graphdir="../graphs/"
print("Ploting the density function graph in "+graphdir+"...")
plt.figure()

#plot contours
plt.contour(lon*rad2deg,lat*rad2deg,data, colors='black')
cp = plt.contourf(lon*rad2deg,lat*rad2deg,data,100,cmap='jet')

#label the axis
plt.xlabel('Longitude', fontsize=14)
plt.ylabel('Latitude', fontsize=14)

#plot colobar
cbar = plt.colorbar(cp)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

#show, save and close
plt.savefig(graphdir+'density_function.png',format='png')
plt.show()
plt.close()
print("Done")

