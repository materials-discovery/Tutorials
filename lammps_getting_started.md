

# LAMMPS on Aristotle@UCL

# Introduction 
This is a getting started tutorial/how to for setting up and running LAMMPS [https://www.lammps.org](https://www.lammps.org) on Aristotle at UCL. 

## Disclaimer 
It is not intended to be a full tutorial for LAMMPS it self, but rather an initial document enabling new students to find their path in the information stream, and also allow more experienced users a place to store links and share tips about LAMMPS, specifically designed and targeting the Student Master Projects.

All Logos and external content is copyrighted with respective owners, all links, are provided for information purposes only as is, if you spot something wrong or inappropriate please contact A.H.

## What is LAMMPS
* LAMMPS official website:  https://www.lammps.org](https://www.lammps.org). 

![https://www.lammps.org](tutorials/Tutorials/figures/lammps_logo.png)

LAMMPS stands for Large-scale Atomic/Molecular Massively Parallel Simulator.

> _LAMMPS is a classical molecular dynamics simulation code with a focus on materials modeling. It was designed to run efficiently on parallel computers. It was developed originally at Sandia National Laboratories, a US Department of Energy facility. The majority of funding for LAMMPS has come from the US Department of Energy (DOE). LAMMPS is an open-source code, distributed freely under the terms of the GNU Public License Version 2 (GPLv2)._ 

## Where can I find more information about LAMMPS
The best source of information is the [LAMMPS home page](https://www.lammps.org) and the [LAMMPS manual](https://docs.lammps.org/Manual.html). LAMMPS is also found on the official [github repository](https://github.com/lammps/lammps).

Also see this tutorials page: [https://lammps.org/tutorials.html](https://lammps.org/tutorials.html)

## Prerequisite 
in order to follow the tutorials and how to's here you should first read the `aristotle_howto.md` in this repo and get aquatinted with using the terminal and linux shell. 

# Setting up LAMMPS on Aristotle 

Please run the following commands in the terminal (copy and paste after the `-bash-4.2$ ` or `%` prompt (tip you can change the prompt look, [check this out](https://wiki.archlinux.org/title/Bash/Prompt_customization) or simpley execute this in a terminal `export PS1="$PWD "` and if you like it, just added to the ~./bashrc file in your home folder) 

if using Aristotle, run the following to set some program paths and configurations:  

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

* **Tip** 
	* cut and paste the above into the file ~/.bashrc in your home folder on Aristotle, then you never need to redo this again (hint, use e.g., `nano` and copy paste the above).
	* Alternatively, use the `.bashrc` and  `.bash_profile` files in the folder `utils` in this repo, move them to your  home folder in Aristotle and then log out and then log in. You are all set.   

Now you should be able to type in the command line the following 

```shell
(aristotle03 ~) > which lmp_aristotle 
```
If you see something like this, you are all set for simple LAMMPS jobs:

```shell
/home/ccaabaa/apps/lammps/29Sep2021_update2/basic-aristotle/gnu-10.2.0/bin/lmp_aristotle
(aristotle03 ~) > 
```
# Run an example 

## Prepare a LAMMPS work folder  
It is a good practice to create a new folder for each project/tutorial. 

In this case we create on Aristotle (using the terminal or WINSCP) a new folder in the home folder. You may chose any folder name you want, here we choose "lammps_tut_1" : 

if you prefer the Terminal, type in the terminal prompt: 

* `cd ` + press the Enter key to return to the main folder of your account (always `cd` without arguments gets you back to the home folder, its like the home button on a phone). 
* `mkdir lammps_tut_1` + Enter to create this folder (use any name you want) 
* `cd lammps_tut_1`  + Enter to enter this new folder. (from now on we will assume an Enter at the end of each command!)

now if you list the contents of your home folder on Aristotle like so 

* `ls ` 

you will see no output, as this is a new folder ;-) 



## Get some LAMMPS examples
we want to use some LAMMPS examples, the **easiest way to get them directly to your account** in Aristotle is to clone the entire LAMMPS repository, this way not only the examples are available, but also the documentation and the source code. Of course you can later delete any folder you do not want using the Unix `rm ` command. 

We shall use some of the official LAMMPS tutorials. 

Inside a terminal (on Aristotle)  (don't forget to press Enter at the end of each line to execute it)

* `cd ` + Enter, to got the home folder, we will clone to the home folder, but you can of course clone to any other folder of your choice.  
* `git clone  https://github.com/lammps/lammps.git` 
  
Now you should have a folder named `lammps` in your home folder. This is the entire lammps source code. 

this will take a while, once the entire repository is cloned, you can browse it using `cd lammps
` and `cd ..` etc.   

### Run the simple Melt Example: 

follow these steps in the Terminal on Aristotle: 

* `cd ` + Enter, to got the home folder
* `cd lammps_tut_1`
* `cp -r ~/lammps/examples/melt/ .` --> this will copy the entire `melt` example folder from the LAMMPS examples to your `lammps_tut_1` folder
* `cd melt`
* `lmp_aristotle -i in.melt`  --> this will actually run LAMMPS interactively, i.e, it will run in the foreground in the terminal, the results will be printed directly to the screen. Here is how the first few lines look like: 

```sh
(aristotle03 melt) > lmp_aristotle -in in.melt > out.melt
[aristotle03.kathleen.ucl.ac.uk:62204] mca_base_component_repository_open: unable to open mca_mtl_psm: libpsm_infinipath.so.1: cannot open shared object file: No such file or directory (ignored)
(aristotle03 melt) > lmp_aristotle -in in.melt
[aristotle03.kathleen.ucl.ac.uk:62221] mca_base_component_repository_open: unable to open mca_mtl_psm: libpsm_infinipath.so.1: cannot open shared object file: No such file or directory (ignored)
LAMMPS (29 Sep 2021 - Update 2)
Lattice spacing in x,y,z = 1.6795962 1.6795962 1.6795962
Created orthogonal box = (0.0000000 0.0000000 0.0000000) to (16.795962 16.795962 16.795962)
  1 by 1 by 1 MPI processor grid
Created 4000 atoms
  using lattice units in orthogonal box = (0.0000000 0.0000000 0.0000000) to (16.795962 16.795962 16.795962)
  create_atoms CPU = 0.001 seconds
Neighbor list info ...
.....
.....

```
**Note there is a warning related to mismatch of libraries, this is due to cross compilation for other machines, it can be ignored safely**

* While it is fine to run short tasks interactively, later, for your projects, you would need to run much longer jobs and if the terminal is closed while the job is still running, or you turn off or put your machine to sleep the job (the lammps application) will stop and crash. Hence instead we now run in the background: 
* run LAMMPS in the background with output stored in lammps.out: 



```sh 
lmp_aristotle -in in.file > out.melt & 
```

the `>` is a redirection, it instructs the shell (the UNIX/Linux CLI programme) to redirect all output that normally is printed on the screen (i.e., terminal) to a file, in our case its called out.melt. The `&` means that the job is going to run in the background, so that even when you logout, close the terminal, or just shut down your machine, the job will continue until completion. 

To check if a job is still running use the `jobs` and the `ps` (process command). 
`jobs` --> shows the running jobs running, the output looks like so: 

```sh
(aristotle03 melt) > jobs
[1]-  Running                 lmp_aristotle -in in.melt > out.melt &
[2]+  Running                 lmp_aristotle -in in2.melt > out2.melt & 
``` 
Note: you can run multiple jobs, as above. 

to bring one of the jobs back to the foreground, type

`fg %#` --> where `#` should be replaced by the job number, so if we want to bring forward the second job we do: 

`fg %2` 

and to send it back to the background, just type

`CTRL + Z ` + `bg` + `Enter'

Another option is to use the 'ps' command, but we will leave it for another time. 

**Congratulations! you run the first lammps simulation**

### Inspect the LAMMPS Input file
Lets inspect the lammps input file `in.melt`. Since this file is already on Aristotle we need some tools either on the machine itself, or we need to copy it to our desktop/laptop and look at it. The same of course holds for the out.melt file, and in general to any file on the remote server (which is Aristotle in our case almost always).

One option is to use [scp or winscp]((https://www.rc.ucl.ac.uk/docs/howto/#windows-winscp)) and copy the file to your desktop/laptop, view, edit, change etc., using your windows tool which you are familiar with (you can use any tool which edits text, not wordpad or MS word as they will alter the format and make it void! 

It is recommended to use a tool like [vscode](https://code.visualstudio.com) which is an open source solution from Microsoft, and can even do remote editing in addition to local files, see [https://code.visualstudio.com/docs/remote/ssh](https://code.visualstudio.com/docs/remote/ssh) for how to set it up for remote editing.


In this tutorial we will edit the file directly on Aristotle i.e., without copying it over. For this we will use the simple text editor program called [nano](https://www.nano-editor.org) (alternatively, if you are familiar with `vim`, or are feeling adventures ;-), you may try to learn [vim](https://opensource.com/article/19/3/getting-started-vim) instead, which is already active by default). `nano` has to be activated first: 

* `module load nano`

now run 

* `nano in.melt`  and use the arrow keys to move around. Identify all the various parts of the input file, and consult the lammps manual for help! 
  

  Another Tip: instead of scp/winscp and a CMD terminal, you can use the more advanced, but easy to use  [MobaXterm](https://www.rc.ucl.ac.uk/docs/howto/#mobaxterm)) program!

# Plotting the thermal properties 

## Option 1: creating plots without copying data to your local machine

this requires some advanced but not too difficult knowledge of the Unix/Linux X11 (X Server) which is a standard for network graphical user interfaces and applications widely use in the scientific community. 

For this, if you are using a mac, you need to get [XQuartz](https://www.xquartz.org/) either by installing directly or using the [homebrew](https://brew.sh) system. If you are using a windows machine, use [Putty](https://www.putty.org/) with Exceed (similar to XQuartz...) . The latter is available in the UCL PC @Anywhere desktop. (Hint if you have a MAC or IPAD or even a windows machine etc, you can still use the standard UCL Desktop by following the instructions here [# Desktop@UCL Anywhere](https://www.ucl.ac.uk/isd/services/computers/remote-access/desktopucl-anywhere) )

check also (for windows) [### Using PuTTY with Exceed](https://www.rc.ucl.ac.uk/docs/Supplementary/X-Forwarding/#using-putty-with-exceed)
 or see summary here [putty wth exceed summary](putty_exeed_ucl.md)
for mac run XQuartz and simply open a terminal and run `ssh -X aritotle.rc.ucl.ac.uk -l <yourusername>` and thats it.

Once you establish an ssh connection to Aristotle using putty+Exceed on windows or terminal + XQuartz on MAC with the X forwarding option (-X on mac, and check the above putty + xceed page for windows) you will be able to run graphical programs remotely. see (UCL Research Computing Documentation)[https://www.rc.ucl.ac.uk/docs/howto/]



so now we can run gnuplot directly on Aristotle: 




## Option 1: creating plots with copying data to your local machine

# Visualisation 


