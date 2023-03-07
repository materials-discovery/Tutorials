#!/bin/bash
#extract thermo data. run this file like so in teh same folder where log.lammps is: 
# ./extract_thermo.sh 
# Dont forget to press enter after you write the command. It will create a new file called thermo.dat which has the thermodynamic data

T=$(cat log.lammps | grep "\s*fix\s*1\s*all\s*nvt" | awk '{print $7}')
n_run=$(cat log.lammps | grep "run" | awk '{print $2}')
cat log.lammps| sed -n '/Step/, /Loop time/p' | awk ' NR == 1 { print "#" $0 }  !/[Step,Loop]/ {print  }' > thermo.dat
