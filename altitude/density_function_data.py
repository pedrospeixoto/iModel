#! /usr/bin/env python3
#---------------------------------
#   Luan Santos (luan.santos@usp.br)
#   September 2021
#----------------------------------

#----------------------------------
# This script generates the density function data
# defined in a lat-lon grid
# Output: densf_table.dat (in altitude directory)
# To use this data in iModel grid generators, you must
# set the following options in mesh.par file:
#
# !Position: eqs, pol, ran, ref, readref, readref_andes
# readref
# !Optimization of mesh: nopt, sprg, scvt, salt, hr95
# scvt
#----------------------------------
import numpy as np
import matplotlib.pyplot as plt

def density_function(lat,lon):
   gammas = 5.0
   
   #radius (in km) of high resolution area
   maxdist = 100.0
   
   #distance (in km) of transition belt
   epsilons = 4000.0
   
   #center of refined region is (0,0)
   #distance to center (Haversine Formula)
   coslat = np.cos(lat)
   radiuse = 6371.0
   dists = radiuse*2.0*np.arcsin(np.sqrt( np.sin(lat/2.0)**2 + coslat*np.sin(lon/2.0)**2  ))
   
   #distance to center metric
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
data = np.zeros((nlat+1,nlon+1))
print("Done!")

#compute the density function in the lat-lon grid
data = density_function(lat,lon)

#plot the data
plt.figure()

#plot contours
#converts degrees to radians
rad2deg = 1.0/deg2rad
plt.contour(lon*rad2deg,lat*rad2deg,data, colors='black')
cp = plt.contourf(lon*rad2deg,lat*rad2deg,data,100,cmap='jet')

#label the axis
plt.xlabel('Longitude', fontsize=16)
plt.ylabel('Latitude', fontsize=16)

#plot colobar
cbar = plt.colorbar(cp)
cbar.formatter.set_powerlimits((0, 0))
cbar.update_ticks()

#show, save and close
plt.savefig('data.png',format='png')
plt.show()
plt.close()

#saves the data in densf_table
print("Saving the data in a file...")
with open('densf_table.dat', 'w') as f:
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
print("Done!")
