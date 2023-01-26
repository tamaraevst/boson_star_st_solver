#!/bin/bash

read -p "Folder with data:" folder 
read -p "Min value for central gravitational field:" phigravmin 
read -p "Max value for central gravitational field:" phigravmax
read -p "Name of your parameter file (pls provide respective path if not in current directory):" paramfile

# folder=beta_12

for i in $(seq $phigravmin 0.01 $phigravmax)
do
   sed -i.bak '25d' $paramfile

   echo "phigrav0           ${i}"  >> $paramfile

   ./../single_condition $paramfile   

   mkdir $folder

   mv phigrav.dat $folder/phigrav_${i}.dat
   mv X.dat $folder/X_${i}.dat
   mv eta.dat $folder/eta_${i}.dat
   mv m.dat $folder/m_${i}.dat
   mv Phi.dat $folder/Phi_${i}.dat
   mv A.dat $folder/A_${i}.dat
   cp $paramfile $folder/.

done

python3 recursive_plot.py -i $folder

