#!/bin/bash
#TC2 SCVT experiments

# Script for Yonggang Yu to reproce results for SW TC2 experiments
#!/bin/bash

source inc_vars.sh 

# Script for convergence analysis on HR95 grids
#------------------------------------------

#Setup enviroment 
make

#Adjust threading
export OMP_NUM_THREADS=4

#Download grids
./sh/getgridsfromserver.sh


echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cp par/mesh_Yonggang.par par/mesh.par
cat par/mesh.par


awk '{ if ( NR == 5 ) { print "read";} else {print $0;} }'  par/mesh.par > par/mesh2.par
awk '{ if ( NR == 9 ) { print "nopt";} else {print $0;} }'  par/mesh2.par > par/mesh.par
#cp par/mesh2.par par/mesh.par

cp par/swm_Yonggang.par par/swm.par
cat par/swm.par



awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_1.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_2.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_3.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_4.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_5.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_6.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_7.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_8.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel

awk '{ if ( NR == 17 ) { print "icos_pol_scvt_h1_9.xyz";} else {print $0;} }'  par/mesh.par > par/mesh2.par
cp par/mesh2.par par/mesh.par
./imodel
