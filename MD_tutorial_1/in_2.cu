#lammps Cu melt
# initialise temperature distribution


units metal
boundary p p p 
atom_style atomic 


# lets input the lattice
read_data 4105040.lmp

# let us define the interatomic potential

pair_style adp 
# map the type 1 atom in lammps input file to the Cu in the potential file, see https://docs.lammps.org/pair_adp.html
pair_coeff * * CuTa_LJ15_2014.adp Cu

neighbor	0.3 bin
neigh_modify	delay 10

fix		1 all nve

velociy all create 1200 2134983 dist gaussian 

# https://docs.lammps.org/thermo.html
thermo		1

timestep	0.001

# https://docs.lammps.org/dump.html
dump		1 all atom 50 dump_cu.adp

#dump		2 all image 10 image.*.jpg element element &
#		axes yes 0.8 0.02 view 60 -30
#dump_modify	2 pad 3 element Si C

#dump		3 all movie 10 movie.mpg element element &
#		axes yes 0.8 0.02 view 60 -30
#dump_modify	3 pad 3 element Si C

run		1000






