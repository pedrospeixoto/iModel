#!/bin/bash

# Script for multiple deformational tests
#------------------------------------------

source inc_vars.sh 

echo "            Basic Parameters     "
echo "---------------------------------"
echo " "
cat deform.par

echo " "
echo "            Test 1 - Case 3 - div velocity - Ini 1-Cosbell  "
echo "------------------------------------------------------------"
echo " "
#Test case
awk '{ if ( NR == 3 ) { print 3;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 1;} else {print $0;} }'  deform2.par > deform.par

cat deform.par

#./driver
# runtests
#./runntimes.sh

#Creat maps
cd ../gmt
./anim.sh ../data/deform_case3_ini1_traj1_nt60_phi_t icos_eqs_nopt_5 
./anim.sh ../data/deform_case3_ini1_traj1_nt60_rho_t icos_eqs_nopt_5
./anim.sh ../data/deform_case3_ini1_traj1_nt60_phi_rho_t icos_eqs_nopt_5
cd ../src

echo " "
echo "            Test 2 - Case 4 - Rotation - Ini 1-Cosbell        "
echo "--------------------------------------------------------------"
echo " "
#Test case
awk '{ if ( NR == 3 ) { print 4;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 1;} else {print $0;} }'  deform2.par > deform.par

cat deform.par
#./driver
# runtests
./runntimes.sh

#Create maps
cd ../gmt
./anim.sh ../data/deform_case4_ini1_traj1_nt60_phi_t icos_eqs_nopt_5 
./anim.sh ../data/deform_case4_ini1_traj1_nt60_rho_t icos_eqs_nopt_5
./anim.sh ../data/deform_case4_ini1_traj1_nt60_phi_rho_t icos_eqs_nopt_5
cd ../src


echo " "
echo "            Test 3- Case 4 rotat - Slot cylinder         "
echo "--------------------------------------------------------- "
echo " "
#Test case
awk '{ if ( NR == 3 ) { print 4;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 3;} else {print $0;} }'  deform2.par > deform.par
cat deform.par

# runtests
./runntimes.sh
#Creat maps
cd ../gmt
./anim.sh ../data/deform_case4_ini3_traj1_nt60_phi_t icos_eqs_nopt_5 
./anim.sh ../data/deform_case4_ini3_traj1_nt60_rho_t icos_eqs_nopt_5
./anim.sh ../data/deform_case4_ini3_traj1_nt60_phi_rho_t icos_eqs_nopt_5
cd ../src


echo " "
echo "            Test 4 - Nondiv case 1 - Ini Gauss-2     "
echo "-----------------------------------------------------"
echo " "
#Test case
awk '{ if ( NR == 3 ) { print 1;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 2;} else {print $0;} }'  deform2.par > deform.par
cat deform.par

# runtests
./runntimes.sh

#Creat maps
cd ../gmt
./anim.sh ../data/deform_case1_ini2_traj1_nt60_phi_t icos_eqs_nopt_5 
./anim.sh ../data/deform_case1_ini2_traj1_nt60_rho_t icos_eqs_nopt_5
./anim.sh ../data/deform_case1_ini2_traj1_nt60_phi_rho_t icos_eqs_nopt_5
cd ../src


echo " "
echo "            Test 5 - Non div case 2 - ini gauss 2  "
echo "---------------------------------------------------"
echo " "
#Test case
awk '{ if ( NR == 3 ) { print 2;} else {print $0;} }'  deform.par > deform2.par
#Initial Condition Field
awk '{ if ( NR == 5 ) { print 2;} else {print $0;} }'  deform2.par > deform.par
cat deform.par

# runtests
./runntimes.sh

#Creat maps
cd ../gmt
./anim.sh ../data/deform_case2_ini2_traj1_nt60_phi_t icos_eqs_nopt_5 
./anim.sh ../data/deform_case2_ini2_traj1_nt60_rho_t icos_eqs_nopt_5
./anim.sh ../data/deform_case2_ini2_traj1_nt60_phi_rho_t icos_eqs_nopt_5
cd ../src
