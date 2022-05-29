# LAMMPS on Aristotle@UCL

# Introduction 
This is a tutorial/howto for setting up and running LAMMPS [https://www.lammps.org](https://www.lammps.org) on Aristotle at UCL. 

## Disclaimer 
It is not intended to be a full tutorial for LAMMPS it self, but rather an initial document enabling new students to find their path in the information stream, and also allow more experienced users a place to store links and share tips about LAMMPS, specifically designed and targeting the Student Master Projects.

All Logos and external content is copyrighted with respective owners, all links, are provided for information purposes only as is, if you spot something wrong or inappropriate please contact A.H.

## What is LAMMPS
* LAMMPS official website:  [https://www.lammps.org](https://www.lammps.org). 

![https://www.lammps.org](./figures/lammps_logo.png)

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel Simulator.

> _LAMMPS is a classical molecular dynamics simulation code with a focus on materials modeling. It was designed to run efficiently on parallel computers. It was developed originally at Sandia National Laboratories, a US Department of Energy facility. The majority of funding for LAMMPS has come from the US Department of Energy (DOE). LAMMPS is an open-source code, distributed freely under the terms of the GNU Public License Version 2 (GPLv2)._ 

## Where can I find more information about LAMMPS
The best source of information is the [LAMMPS home page](https://www.lammps.org) and the [LAMMPS manual](https://docs.lammps.org/Manual.html). LAMMPS is also found on the official [github repository](https://github.com/lammps/lammps).

## Prerequisite 
in order to follow the tutorials and howto's here you should first read the `aristotle_howto.md` in this repo and get aquatinted with using the terminal and linux shell. 

# Setting up Lammps on Aristotle 

Please run the following commands in the terminal (copy and paste after the

`-bash-4.2$ ` 

if using Aristotle). 

```bash
module load default-modules
module -f unload compilers mpi gcc-libs
module load beta-modules
module use --append /home/ccaabaa/lib/modulefiles/development /home/ccaabaa/lib/modulefiles/bundles /home/ccaabaa/lib/modulefiles/beta
module load gcc-libs/10.2.0
module load compilers/gnu/10.2.0
module load numactl/2.0.12
module load binutils/2.36.1/gnu-10.2.0
module load ucx/1.9.0/gnu-10.2.0
module load mpi/openmpi/4.0.5/gnu-10.2.0
module load python3/3.9-gnu-10.2.0-aristotle


module load lammps/29sep21up2/basic/gnu-10.2.0-aristotle
```

* **Tip** cut and paste the above into the file .bashrc in your home folder on Aristotle, then you never need to redo this again (hint, use `nano`). A ready to use .bashrc is in the folder utils in this repo, move it to your home folder in Aristotle and then log out and then log in. 


### Temporary solution for getting LAMMPS on ARISTOTLE
Until the module procedure is fixe, you should clone the lammps-tutorial repository from https://github.com/materials-discovery/Lammps_tutorials.git like so: 

* Assuming you are already logged into Aristotle then type the following:

```git clone https://github.com/materials-discovery/Lammps_tutorials.git```


* this will clone this repository to your local folder names Lammps_tutorials. 

* now enter into the lammps_tutorials folder: 

```cd Lammps_tutorials```

* there is a file called `lmp_aristotle.tgz`, which is a tar archive, lets extract it: 

``` tar -xzvf lmp_aristotle.tgz```

Now you will have an executable file called `lmp_aristotle` in the `Lammps_tutorials folder`. However, if you try to run it by writing the command 

`lmp_aristotle' then pressing Enter, you get an error indicating some missing linraries, this is because you need still to load the compiler into your environment, so do the following: 

`module load gcc-libs/4.9.2`

then you can type the `lmp_aristotle` and it will work, it will show something like below: 

![](lmp_check.png)

### Prepare a folder for working on the tutorials 
It is a good practice to create a new folder for each project/tutorial. In this case reate on Aristotle (using the terminal) a new folder in the home or main folder named "lammps_tutorial_1" like so: 

type: 

* `cd ` + Enter to return to the main folder of your account. 
* `mkdir lammps_tutorial_1` + Enter to create 
* `cd lammps_tutorial_1`  + Enter to enter this new folder. (from now on we will assume an Enter at the end of each command!)

now if you list the contents of your home folder on Aristotle like so 



### Clone the LAMMPS GitHub Repository
we want to use some LAMMPS examples, the easiest way to get them directly to your account in Aristotle is to clone the entire LAMMPS repository, this way not only the examples are available, but also the documentation and the source code. Of course you can later delete any folder you do not want using the `rm ` command. 

We shall use some of the official LAMMPS tutorials. 

Do in the terminal within the folder (press Enter at the end of each line to execute it)

* `cd ` + Enter, to got the home folder
* `git clone  https://github.com/lammps/lammps.git` 

this will take a while, once the entire repository is cloned, you can browse it using `cd <folder>
` and `cd ..` etc. 

### Run the simple Melt Example: 

do the following
* `cd ` + Enter, to got the home folder
* `cd lammps_tutorial_1`
* `cp -r ~/lammps/examples/melt/ .` --> this will copy the entire melt folder from the lammps examples to your lammps_tutorial_1 folder
* `cd melt`
* ~/Lammps_tutorials/lmp_aristotle -i in.melt`

**Congratulations! you run the first lammps simulation**

### Inspect the LAMMPS Input file
Lets inspect the lammps input file `in.melt`. 

One option is to use [scp or winscp]((https://www.rc.ucl.ac.uk/docs/howto/#windows-winscp)) and copy the file to your desktop/laptop, view, edit, change etc., then use scp or winscp but in this tutorial we will edit the file directly on Aristotle using the simple text editor program called [nano](https://www.nano-editor.org) (if you are familiar with it, or are feeling adventures ;-), you may try to learn [vim](https://opensource.com/article/19/3/getting-started-vim)). 

First we activate nano (we need to do this every time we log in to Aristotle):

* `modeule load nano`

now run 

* `nano in.melt`  and use the arrow keys to move arround. Identify all the various parts of the input file, and consult the lammps manual for help! 
  

  Tip: instead of scp/winscp and a CMD terminal, you can use the more advanced, but easy to use  [MobaXterm](https://www.rc.ucl.ac.uk/docs/howto/#mobaxterm)) program!



## Visualisation of the results 
In order to visualise the results we need to use our laptop/desktop, as remote visualisation may be inefficient. 
* on your laptop/Desktop, if you are windows user, 
