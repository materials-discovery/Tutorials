# Introduction 
This is a tutorial for setting up and running lammps on Aristotle 
## What is Aristotle
Aristotle is an interactive, Linux-based compute service for teaching, running on three each with 64 gigabytes of RAM and 16 cores. 
The nodes run the Red Hat Enterprise Linux operating system (RHEL 7)
see [ISD](
 https://www.rc.ucl.ac.uk/docs/Other_Services/Aristotle/) for more information.

## How to connect to Aristotle 
See the official documentation on [UCL-IDS](https://www.rc.ucl.ac.uk/docs/howto/) HOW TO. 

# Step1: Setting up Lammps on Aristotle 
Lammps is one of the supported and already available packages but it needs to be activated. once you are logged in execute the following commands: 

 `module load gcc-libs/4.9.2`


`module load compilers/intel/2018/update3`

`module load mpi/intel/2018/update3`

`module load lammps/7Aug19/userintel/intel-2018`


Note: run each command separately by copy and pasting the entire line and hitting enter afterwards. If there are issues, please submit an issue or email me directly. 

Each of these commands will enable various parts of the environment needed to run Lammps (e.g., defining which compilers, libraries and executable to load into your path)

Once these commands have run, you should have the executable `lmp_default` in your path and you are ready to move to the next section. 

# Step 2: downloading the LAMMPS tutorials 
We shall use some of the official LAMMPS tutorials. First we need to download them either 

1. from the [official repository](https://github.com/lammps/lammps/tree/master) using a web browser to your local drive and then using scp to move them to Aristotle (see [instructions here](https://www.rc.ucl.ac.uk/docs/howto/#scp)) or 
2. use wget 
3. use git clone
4. simple we copy them from the install folder on Aristotle which for the time of writing of this how to is in '/shared/ucl/apps/lammps/7Aug2019/USERINTEL/intel-2018/lammps/examples'





We need the following tutorials: 

1. melt:     rapid melt of 3d LJ system
2. min:      energy minimization of 2d LJ melt
3. obstacle: flow around two voids in a 2d channel
4. flow:     Couette and Poiseuille flow in a 2d channel


WIP