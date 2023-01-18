# Prerequisites: 
- [[MD model creation]]
- [[lammps_getting_started]]
Make sure you read and perform the [MD model creation](../MD_Model_Creation/MD model creation.md) tutorial as this one picks up from it. 
We assume also you have lammps already cloned in your Aristotle account under the folder `/home/<username>/lammps` (<> mean replace whats inside with username etc).

On **Aristotle** (this tutorial assumes you are using it) do the following: 

1. make a new working folder called for example:
`mkdir lammps_Cu_lattpar`
and change folder into it: 
`cd lammps_Cu_lattpar`
now when you  run the command `pwd` you should see that you are in this folder. 
2. let us copy the Cu we made in  [MD model creation](../MD_Model_Creation/MD model creation.md)  to this folder, assuming it is on the same level: 
`cp ../MD\ model\ creation/4105040.lmp .`
note that the `..` refers to the parent folder, and the `\ ` that is backslash and a space is the space character in unix, and the last ` . ` is a shortcut for the current folder, i.e., we do a copy (cp) of the file `4105040.lmp` from the path: `../MD\ model\ creation/4105040` to the current folder. 

3. Now we could use our favorite method (`nano`, `vi` or `winscp`,  `mobaxterm` etc) to create an input file for LAMMPS, but it is tedious to create it from scratch, so we use this one as a reference, from the LAMMPS example folder (see  [[lammps_getting_started]]). This is a good way to start a LAMMPS input file in general. 
	1. tip check the lammps/examples/Readme file to check what examples are available. 
	2. we will use as a simple example,  `lammps/examples/meam/in.meam` which uses  a modified embedded atom (MAEM) potential, which is as the name suggests, similar to the EAM potential we learned about in the lectures.
	
```js
# Test of MEAM potential for SiC system
units metal
boundary p p p
atom_style atomic
read_data data.meam    # L1 will be modified
pair_style meam        # L2 will be modified
pair_coeff * * library.meam Si C SiC.meam Si C # L3 will be modified
neighbor 0.3 bin 
neigh_modify delay 10
fix 1 all nve
thermo 10
timestep 0.001

#dump 1 all atom 50 dump.meam
#dump 2 all image 10 image.*.jpg element element &
# axes yes 0.8 0.02 view 60 -30
#dump_modify 2 pad 3 element Si C
#dump 3 all movie 10 movie.mpg element element &
# axes yes 0.8 0.02 view 60 -30
#dump_modify 3 pad 3 element Si C
run 100

``` 
copy paste this potential into a text file in your working folder `lammps_Cu_lattpar` using e.g., `nano` as we learned in the tutorial
`% nano in_1.cu` 
and use Ctrl+C and Ctrl+V to paste it, then save the file and exit. 
if you do an ls you should see: 
```
% ls
in_1.cu
```


This is your first input file for LAMMPS (unless of course you already have done it in the past!)

Below we will have to change the lines marked with L1, L2, and L3 to reflect the new potential (you can also use the EAM potential as is if you like, but its for a different material)

4. Let us find and download an interatomic potential from the NIST repository [NIST Interatomic Potential Repository](https://www.ctcms.nist.gov/potentials/)
	1. we choose the CuTa potential (a newer version of the one we have seen in the first MD lecture). 
	2. It contains Cu-Cu, Cu-Ta and Ta-Ta interactions, we will only use now the Cu-Cu interactions. Note, that one should always test the properties for known properties, especially if we are using a binary alloy potential for one element, as the fitting process may have preferred the CU-Ta interaction over the Cu-Cu or Ta-Ta as one in these systems is more focused on the interfaces or alloys, but this potential is known to have good representation of each of the phases. But we will test it now for simple property, the lattice parameter. 
5. change the inout file so that 
	1. L1 is replaced with `read_data 4105040.lmp`
		1. note that `4105040.lmp` is the file we created in the previous tutorial. 
	2. L2 is replaced with: `pair_style adp ` 
	3. finally L3 is replace by: `pair_coeff * * CuTa_LJ15_2014.adp Cu`
4. your Input file should like something like this:  (note we added some comments, and the order may be different)
```js
#lammps Cu melt


units metal
boundary p p p 
atom_style atomic 


# lets input the lattice
read_data 4105040.lmp

# let us define the interatomic potential

pair_style adp 
# map the type 1 atom in lammps input file to the Cu in the potential file, see https://docs.lammps.org/pair_adp.html
pair_coeff * * CuTa_LJ15_2014.adp Cu

neighbor        0.3 bin
neigh_modify    delay 10

fix             1 all nve
# https://docs.lammps.org/thermo.html
thermo          1

timestep        0.001

# https://docs.lammps.org/dump.html
dump            1 all atom 50 dump_cu.adp

#dump           2 all image 10 image.*.jpg element element &
#               axes yes 0.8 0.02 view 60 -30
#dump_modify    2 pad 3 element Si C

#dump           3 all movie 10 movie.mpg element element &
#               axes yes 0.8 0.02 view 60 -30
#dump_modify    3 pad 3 element Si C

run             1000
```

5. now you are ready to run it:
`% lmp < in_1.cu > out_1.cu`

the output file `out_1.cu` should look like this
```js
LAMMPS (24 Mar 2022)
Reading data file ...
  orthogonal box = (0 0 0) to (17.90955 17.90955 17.90955)
  1 by 1 by 1 MPI processor grid
  reading atoms ...
  365 atoms
  read_data CPU = 0.091 seconds
Neighbor list info ...
  update every 1 steps, delay 10 steps, check yes
  max neighbors/atom: 2000, page size: 100000
  master list distance cutoff = 6.450959
  ghost atom cutoff = 6.450959
  binsize = 3.2254795, bins = 6 6 6
  1 neighbor lists, perpetual/occasional/extra = 1 0 0
  (1) pair adp, perpetual
      attributes: half, newton on
      pair build: half/bin/atomonly/newton
      stencil: half/bin/3d
      bin: standard
Setting up Verlet run ...
  Unit style    : metal
  Current step  : 0
  Time step     : 0.001
Per MPI rank memory allocation (min/avg/max) = 6.776 | 6.776 | 6.776 Mbytes
   Step          Temp          E_pair         E_mol          TotEng         Press     
         0   0             -1169.9166      0             -1169.9166     -22822.181    
        10   2.1887862     -1170.0196      0             -1169.9166     -23197.297    
        20   8.352394      -1170.3097      0             -1169.9167     -24294.29     
        30   17.406901     -1170.7357      0             -1169.9167     -26037.093    
        40   27.907964     -1171.2299      0             -1169.9168     -28325.614    
        50   38.453263     -1171.726       0             -1169.9168     -31061.93     
        60   48.053041     -1172.1777      0             -1169.9168     -34173.181    
        70   56.347388     -1172.568       0             -1169.9168     -37624.526    
        80   63.539027     -1172.9063      0             -1169.9168     -41400.899    
        90   70.01953      -1173.2112      0             -1169.9168     -45488.48     
       100   75.919789     -1173.4889      0             -1169.9168     -49844.434    
Loop time of 0.0517902 on 1 procs for 100 steps with 365 atoms

Performance: 166.827 ns/day, 0.144 hours/ns, 1930.867 timesteps/s
98.7% CPU use with 1 MPI tasks x no OpenMP threads

MPI task timing breakdown:
Section |  min time  |  avg time  |  max time  |%varavg| %total
---------------------------------------------------------------
Pair    | 0.050085   | 0.050085   | 0.050085   |   0.0 | 96.71
Neigh   | 0.00024979 | 0.00024979 | 0.00024979 |   0.0 |  0.48
Comm    | 0.0003455  | 0.0003455  | 0.0003455  |   0.0 |  0.67
Output  | 0.00081796 | 0.00081796 | 0.00081796 |   0.0 |  1.58
Modify  | 0.00019812 | 0.00019812 | 0.00019812 |   0.0 |  0.38
Other   |            | 9.421e-05  |            |       |  0.18

Nlocal:            365 ave         365 max         365 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Nghost:           1336 ave        1336 max        1336 min
Histogram: 1 0 0 0 0 0 0 0 0 0
Neighs:          12007 ave       12007 max       12007 min
Histogram: 1 0 0 0 0 0 0 0 0 0

Total # of neighbors = 12007
Ave neighs/atom = 32.89589
Neighbor list builds = 1
Dangerous builds = 0
Total wall time: 0:00:00




```

6. the above runs only 100 steps, which as we discussed in class, does not lead to significantly relevant statistics, in the class tutorial we shall 
	1. increase the number of steps 
	2. initialise the velocities to a desired temperature and distribution (using the command `velociy all create 1200 2134983 dist gaussian `)
	3. change the ensemble to NVT to bring the system into specific target temperatures 
	4. change the ensemble to NPT to reach a specific temperature and pressure
	5. In all cases we will analyse the output log file and plot the total energy and temperature as a function of time (snapshots) and we will also visualise the system.
	6. relax the system to find the lattice parameter using the `https://docs.lammps.org/fix_box_relax.html` command
		1. Hint we remove the dynamic run and 