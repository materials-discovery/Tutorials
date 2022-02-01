# LAMMPS on Aristotle@UCL

## Introduction 
This is a tutorial for setting up and running Lammps on Aristotle at UCL. 

## What is LAMMPS
![Drag Racing](lammps-logo.png)

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel Simulator.

LAMMPS is a classical molecular dynamics simulation code with a focus on materials modeling. It was designed to run efficiently on parallel computers. It was developed originally at Sandia National Laboratories, a US Department of Energy facility. The majority of funding for LAMMPS has come from the US Department of Energy (DOE). LAMMPS is an open-source code, distributed freely under the terms of the GNU Public License Version 2 (GPLv2).

## Where can I find more information about LAMMPS
The best source of information is the [LAMMPS home page](https://www.lammps.org) and the [LAMMPS manual](https://docs.lammps.org/Manual.html). LAMMPS is also found on a [github repository](https://github.com/lammps/lammps) 

## What is Aristotle
Aristotle is an interactive, Linux-based compute service for teaching, running on three each with 64 gigabytes of RAM and 16 cores. 
The nodes run the Red Hat Enterprise Linux operating system (RHEL 7)
see [ISD](
 https://www.rc.ucl.ac.uk/docs/Other_Services/Aristotle/) for more information.

## How to connect to Aristotle 
See the official documentation on [UCL-IDS](https://www.rc.ucl.ac.uk/docs/howto/) HOW TO. 

## You need to use a Terminal for this tutorial
For this tutorial you need a Terminal program on your desktop and connect to a shell (a remote terminal) in Aristotle. ![](Linux_command-line._Bash._GNOME_Terminal._screenshot.png).

In windows 10/11 you can either use [putty](https://www.rc.ucl.ac.uk/docs/howto/#windows-putty) or open a cmd or powershell (or WSL2) while for Mac/Linux a local terminal then use ssh to connect to Aristotle. For more information and help about the linux/maac (in fact, more generally UNIX) Terminal see [A guide to the Linux terminal for beginners](https://opensource.com/article/21/8/linux-terminal) and a nice intro to the [main commands you need in a terminal](https://www.redhat.com/sysadmin/10-commands-terminal).

# Step1: Setting up Lammps on Aristotle 
**NOTE: Please do not yet try this as there is a bug on Aristotle and our Research Computing devision is taking care of it, in the mean time, please copy an executable from this repo. **

Lammps is one of the supported and already available packages but it needs to be activated. Once you are logged in execute the following commands: 

 `_> module load gcc-libs/4.9.2`


`_> module load compilers/intel/2018/update3`

`_> module load mpi/intel/2018/update3`

`_> module load lammps/7Aug19/userintel/intel-2018`


Note: run each command separately by copy and pasting the entire line and hitting enter afterwards. If there are issues, please submit an issue or email me directly. 

Each of these commands will enable various parts of the environment needed to run Lammps (e.g., defining which compilers, libraries and executable to load into your path)

Once these commands have run, you should have the executable `lmp_default` in your path and you are ready to move to the next section. 

# Step 2: Prepare a folder for working on the tutorials 
It is a good practice to create a new folder for each project/tutorial. In this case reate on Aristotle (using the terminal) a new folder named "lammps_tutorial" like so: 

`_> mkdir lammps_tutorial`

now if you list the contents of your home folder on Aristotle like so 

`_> ls -lrtha`

you will see a list of folders and files, including the new `lammps_tutorial` folder that you just created. 
![](Screenshot%202022-01-31%20at%2021.34.30.png)


# Step 3: Copy or download the LAMMPS tutorials
We shall use some of the official LAMMPS tutorials. First we need to download them either 

1. from the [official repository](https://github.com/lammps/lammps/tree/master) using a web browser to your local drive and then using scp to move them to Aristotle (see [instructions here](https://www.rc.ucl.ac.uk/docs/howto/#scp)) 
OR 
2. use git clone to a local folder of the entire LAMMPS github repository including source code like so (this is run in a terminal on Aristotle): `_> git clone git@github.com:lammps/lammps.git`
OR 
3. simply, we copy them from the install folder available on Aristotle which for the time of writing of this tutorial is in `/shared/ucl/apps/lammps/7Aug2019/USERINTEL/intel-2018/lammps/examples`

**In this tutorial we follow option 3.** 

## Example 1:  Atomistic [Couette flow](https://en.wikipedia.org/wiki/Couette_flow) and [Poiseuille flow](https://en.wikipedia.org/wiki/Hagen–Poiseuille_equation#Plane_Poiseuille_flow) 

In the terminal, while logged in Aristotle perform the following: 

`_> cp -r /shared/ucl/apps/lammps/7Aug2019/USERINTEL/intel-2018/lammps/examples/flow ./`

This is a copy command (cp) that copies the flow case from the lammps example folder (the `-r` options stands for copying all files/content recursively) to a local flow  folder in your lammps_tutorial on Aristotle. 

Lets inspect the lammps input file `in.flow.couette`. One option is to use scp and copy the file to your desktop, view, edit, change etc., then use scp to copy back (there are other options e.g. using [Windows winscp](https://www.rc.ucl.ac.uk/docs/howto/#windows-winscp) or [MobaXterm](https://www.rc.ucl.ac.uk/docs/howto/#mobaxterm)) but in this tutorial we will edit the file directly on Aristotle using the simple text editor program called [nano](https://www.nano-editor.org) (if you are familiar with it, or are feeling adventures ;-), you may try to learn [vim](https://opensource.com/article/19/3/getting-started-vim)). 

First we activate nano (we need to do this every time we log in to Aristotle):

`_> modeule load nano`

now run 

`_> nano in.flow.couette `

this will open the input file, use the arrow keys to navigate, and the Ctrl+X to exit (say y for saving when promted or n to cancel changes).

In class we will identify the various sections, and will use the [LAMMPS manual](https://docs.lammps.org/Manual.html). 

(Hint you can quickly view a file by the command `_> more <filename>`)

```# 2-d LJ flow simulation

dimension	2
boundary	p s p

atom_style	atomic
neighbor	0.3 bin
neigh_modify	delay 5

# create geometry

lattice		hex 0.7
region		box block 0 20 0 10 -0.25 0.25
create_box	3 box
create_atoms	1 box

mass		1 1.0
mass		2 1.0
mass		3 1.0

# LJ potentials

pair_style	lj/cut 1.12246
pair_coeff	* * 1.0 1.0 1.12246

# define groups

region	     1 block INF INF INF 1.25 INF INF
group	     lower region 1
region	     2 block INF INF 8.75 INF INF INF
group	     upper region 2
group	     boundary union lower upper
group	     flow subtract all boundary

set	     group lower type 2
set	     group upper type 3

# initial velocities

compute	     mobile flow temp
velocity     flow create 1.0 482748 temp mobile
fix	     1 all nve
fix	     2 flow temp/rescale 200 1.0 1.0 0.02 1.0
fix_modify   2 temp mobile

# Couette flow

velocity     lower set 0.0 0.0 0.0
velocity     upper set 3.0 0.0 0.0
fix	     3 boundary setforce 0.0 0.0 0.0
fix	     4 all enforce2d

# Poiseuille flow

#velocity     boundary set 0.0 0.0 0.0
#fix	     3 lower setforce 0.0 0.0 0.0
#fix	     4 upper setforce 0.0 NULL 0.0
#fix	     5 upper aveforce 0.0 -1.0 0.0
#fix	     6 flow addforce 0.5 0.0 0.0
#fix	     7 all enforce2d

# Run

timestep	0.003
thermo		500
thermo_modify	temp mobile

#dump		1 all atom 500 dump.flow

#dump		2 all image 100 image.*.jpg type type &
#		zoom 1.6 adiam 1.2
#dump_modify	2 pad 5

#dump		3 all movie 100 movie.mpg type type &
#		zoom 1.6 adiam 1.2
#dump_modify	3 pad 5

run		10000

```

Let us consider some of the main commands more closely: 

* [lattice command](https://docs.lammps.org/lattice.html)
* [boundary command](https://docs.lammps.org/boundary.html)
* [compute command](https://docs.lammps.org/compute.html)
* [fix command](https://docs.lammps.org/fix.html)

check [LAMMPS manual](https://docs.lammps.org/Manual.html) for others. 

### Run the simulation
`_> 

Note, there seems to be missmatch between the compilers installed and the actual version of lammps, hence we shall now compile lammps from source code. 

`_> mkdir lammps`

`\> cd lammps`

`_> wget --no-check-certificate https://download.lammps.org/tars/lammps-stable.tar.gz `


this will download a file called lammps-stable.tar.gz which is a compressed tar archive, lets uncompress and extract it:

`_> tar -xzvf lammps-stable.tar.gz`
