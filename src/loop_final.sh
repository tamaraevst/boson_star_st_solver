# Script to search the parameter space for BS solution in ST theory of gravity

#!/bin/bash
set -e
set -x

######################################################################
#---------------------------------------------------------------------
######################################################################

# Looped over parameters

min_omega=0.65
max_omega=1.8
N3=31
domega=$(echo "(${max_omega}-${min_omega})/(${N3}-1)"| bc -l)

min_A=0.15788
max_A=0.2
N=21
dA=$(echo "(${max_A}-${min_A})/(${N}-1)"| bc -l)

min_phigrav=-0.1
max_phigrav=0.1
N2=21
dphigrav=$(echo "(${max_phigrav}-(${min_phigrav}))/(${N2}-1)"| bc -l)

######################################################################

# Fixed parameters

nint=128001
rmax=1200

potential="series"
lambda4=-100
lambda6=2500
lambda8=0

thresh=1e-15
minmax="min"
rmatchfac=0.9
mpercentage=99
nzerotarget=0
docheck=0
verbose=0

mphigrav0=0.0
alpha0=0.001
beta0=-10
#omega0=1

######################################################################
#---------------------------------------------------------------------
######################################################################
# Prep work

rm -f zzdata.in
rm -f zzresult.dat


######################################################################

echo "#boson_potential ${boson_potential}" >> zzresult.out
echo "#scalar_potential ${scalar_potential}" >> zzresult.out
echo "#mphigrav0 ${mphigrav0}" >> zzresult.dat
echo "#alpha0 ${alpha0}" >> zzresult.dat
echo "#beta0 ${beta0}" >> zzresult.dat
echo "#A0   phigrav0   omega   star mass   compactness   radius" >> zzresult.dat

for ((i=1;i<$N;i++));
do
	for ((j=1;j<$N2;j++));
	do

		for ((k=1;k<$N3;k++));
		do

					touch zzdata.in
					omega_ARGS=$(echo "$min_omega+($k-1)*$domega"| bc -l)
                	A_ARGS=$(echo "$min_A+($i-1)*$dA"| bc -l)
		            phigrav_ARGS=$(echo "$min_phigrav+($j-1)*$dphigrav"| bc -l)

					echo ""
					echo "START"

					echo "omega  = ${omega_ARGS}"
					echo "min_phigrav = ${dphigrav}"
					echo "A0   = ${A_ARGS}"
					echo "phigrav0   = ${phigrav_ARGS}"
					echo "A0 ${A_ARGS}" >> zzdata.in
					echo "phigrav0 ${phigrav_ARGS}" >> zzdata.in
					echo "mphigrav0 ${mphigrav0}" >> zzdata.in
					echo "alpha0 ${alpha0}" >> zzdata.in
					echo "beta0 ${beta0}" >> zzdata.in
					echo "omega0 ${omega_ARGS}" >> zzdata.in
					echo "nint ${nint}" >> zzdata.in
					echo "rmax ${rmax}" >> zzdata.in
					echo "potential ${potential}" >> zzdata.in
					echo "lambda4 ${lambda4}" >> zzdata.in
					echo "lambda6 ${lambda6}" >> zzdata.in
					echo "lambda8 ${lambda8}" >> zzdata.in
					echo "thresh ${thresh}" >> zzdata.in
					echo "minmax ${minmax}" >> zzdata.in
					echo "rmatchfac ${rmatchfac}" >> zzdata.in
					echo "mpercentage ${mpercentage}" >> zzdata.in
					echo "zero crossings ${nzerotarget}" >> zzdata.in
					echo "docheck ${docheck}" >> zzdata.in
					echo "verbose ${verbose}" >> zzdata.in
					
					#May need to change the executable here appropriately. 
					./single_massless zzdata.in || :
		
                    if test -f "final_model_data_massless.dat"; then
                        sed -n '3p' final_model_data_massless.dat >> zzresult.dat
                    fi
		done
	done
done

rm A.dat Phi.dat Psi0.dat X.dat eta.dat final_model_data.dat joint_profiles.dat m.dat phigrav.dat

mv zzresult.dat zzdone.dat
