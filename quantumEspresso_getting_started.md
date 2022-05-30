# Quantum ESPRESSO on Aristotle@UCL

# Introduction 
This is a getting started tutorial/howto for setting up and running quantum Espresso [https://www.quantum-espresso.org](https://www.quantum-espresso.org) on Aristotle at UCL. 

## Disclaimer 
It is not intended to be a full tutorial for Quantum ESPRESSO it self, but rather an initial document enabling new students to find their path in the information stream, and also allow more experienced users a place to store links and share tips about Quantum ESPRESSO, specifically designed and targeting the Student Master Projects. This is essentially a getting started document. 

All Logos and external content is copyrighted with respective owners, all links, are provided for information purposes only as is, if you spot something wrong or inappropriate please contact A.H.

## What is QE
* Quantum ESPRESSO official website:  [https://www.quantum-espresso.org](https://www.quantum-espresso.org). 

![https://www.quantum-espresso.org](./figures/quantum_ogo_ok.png)

> Quantum ESPRESSO is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.



## Where can I find more information about QE
The best source of information is the [Quantum ESPRESSO home page](https://www.quantum-espresso.org/) and the [documentation](https://www.quantum-espresso.org/documentation/). Quantum ESPRESSO is also found on the official [github repository](https://github.com/QEF/q-e).

Also see the official documentation and tutorials page: [https://www.quantum-espresso.org/tutorials/](https://www.quantum-espresso.org/tutorials/)

## Prerequisite 
in order to follow the tutorials and howto's here you should first read the `aristotle_howto.md` in this repo and get aquatinted with using the terminal and linux shell. 

# Setting up QE on Aristotle 

Please run the following commands in the terminal (copy and paste after the

`-bash-4.2$ ` 

on Aristotle). 

```bash
module load default-modules
module -f unload compilers mpi gcc-libs
module load beta-modules
module use --append /home/ccaabaa/lib/modulefiles/development /home/ccaabaa/lib/modulefiles/bundles /home/ccaabaa/lib/modulefiles/beta
module load gcc-libs/10.2.0
module load compilers/intel/2018/update3
module load mpi/intel/2018/update3
module load xorg-utils/X11R7.7
module load gnuplot/5.0.1
module load quantum-espresso/6.5-impi/thermo_pw-1.2.1/intel-2018
```

* **Tip** 
* 
* cut and paste the above into the file .bashrc in your home folder on Aristotle, then you never need to redo this again (hint, use `nano`). alternatively, use the `.bashrc` and `.bash_profile' files in the folder `utils` in this repo, move them to your  home folder in Aristotle and then log out and then log in.  Note, this will give you also the option to run LAMMPS too, see [LAMMPS getting Started](./lammps_getting_started.md).

Now you should be able to type in the command line the following for the plane wave pseudo potential code: 

```shell
(aristotle03 ~) > which pw.x  
```
If you see something like this, you are all set for simple lammps jobs:

```shell
(aristotle01 ~) > which pw.x
/shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin/pw.x
(aristotle01 ~) > 
```
# Run an example 

## Prepare a QE work folder  
It is a good practice to create a new folder for each project/tutorial. In this case we create on Aristotle (using the terminal) a new folder in the home folder. You may chose any folder name you want, here we choose "qe_tut_1" : 

type in the terminal prompt: 

* `cd ` + press the Enter key to return to the main folder of your account (always `cd` without arguments gets you back to the home folder, its like the home button on a phone). 
* `mkdir qe_tut_1` + Enter to create this folder (use any name you want) 
* `cd qe_tut_1`  + Enter to enter this new folder. (from now on we will assume an Enter at the end of each command!)

now if you list the contents of your home folder on Aristotle like so 

* `ls ` 

you will see no output, as this is a new folder ;-) 



## Get some QE examples
we want to use some QE examples, the easiest way to get them directly to your account in Aristotle is to clone the entire QE repository, this way not only the examples are available, but also the documentation and the source code. Of course you can later delete any folder you do not want using the Unix `rm ` command. 

Inside a terminal (on Aristotle)  (don't forget to press Enter at the end of each line to execute it)

* `cd ` + Enter, to got the home folder, we will clone to the home folder, but you can of course clone to any other folder of your choice.  

We need a specific version of QE. The developers use GitHub tags to mark specific branches and released, so we use a slightly different git clone: 

* `git clone --depth 1 --branch <tag_name> <repo_url>  ` 
* where we replace <tag-name> with the version supported on Aristotle, in our case its `qe-6.5` and the `<repo_url>` is simply:  `https://github.com/QEF/q-e.git` (note we use the HTTPS protocol as we are not members of the dev team, i.e, they can use the git protocol...) so we need to run in Aristotle Terminal: 

```sh
git clone --depth 1 --branch qe-6.5 https://github.com/QEF/q-e.git
```
  
Now you should have a folder named `q-e` in your home folder. This is the entire q-e source code. 

this will take a while, once the entire repository is cloned, you can browse it using `cd q-e
` and `cd ..` etc.   

### Run a simple example: 

follow these steps in the Terminal on Aristotle: 

* `cd ` + Enter, to got the home folder
* `cd qe_tut_1`
* `cp -r ~/q-e/PP/examples/example01.`  --> this will copy the entire `example01` folder from the qe_tur_1 folder


Check the `README` file inside the example01 folder describing the various steps the script is performing. The `run_example` file is a shell script (like .bat files in windows) which prepares input files, executes them, and does some pre and post processing. It is a good practice to open the file and inspect it using e.g. `nano`, or use windSCP to copy and view, see [file_howto.md](aristotle_howto.md) for more information. 

you can also check the top level readme file of the example folder of the Plane wave code found in the path: `~/q-e/PP/examples` which tells us that `example01` shows how to use pw.x and postprocessing codes to make a contour plot in the [110] plane of the charge density for Si, and to plot the band structure of Si.

## Environment variables 
Before starting run_example, we need to set some variables in the file 'environment_variables' so lets open it in Aristotle using `nano` and change the following

1. change the PREFIX line (line 20) to 
  ```
  PREFIX=`cd /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/ ; pwd`
  ```
2. change `TMPDIR` from `TMP_DIR=$PREFIX/tempdir`  to `TMP_DIR=/tmp/<UCL-USERNAME>_tempdir`
Replace <UCL-USERNAME>  of course! 

The resulting part of the file should look like so: 
```sh

PREFIX=`cd /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/ ; pwd`
BIN_DIR=$PREFIX/bin
PSEUDO_DIR=$PREFIX/pseudo
# Beware: everything in $TMP_DIR will be destroyed !
#TMP_DIR=$PREFIX/tempdir
TMP_DIR=/tmp/<ucl-username>_tempdir
```

Now, you can also clone the file from this repo under utils ;-) then you need to copy it instead of the above file and put it in the same place. 

Now, you can go into the example 01 folder and run the example like so (type the ollowing in a terminal) 

```sh 
cd ~/qe_tut_1/example01/
sh run_example
```

If all is ok so far, you should see an output like so: 
```sh
(aristotle01 example01) > sh run_example 
/home/ucl-user/test/q-e/PP/examples/example01 : starting

This example shows how to use pw.x and postprocessing codes to make a
contour plot in the [110] plane of the charge density for Si, and to
plot the band structure of Si.

  executables directory: /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin
  pseudo directory:      /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/pseudo
  temporary directory:   /tmp/ucl-user_tempdir
  checking that needed directories and files exist... done

  running pw.x as:       mpirun -np 4 /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin/pw.x  -nk 1 -nd 1 -nb 1 -nt 1 
  running pp.x as:       mpirun -np 4 /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin/pp.x  -nk 1 -nd 1 -nb 1 -nt 1 
  running plotrho.x as:  /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin/plotrho.x
  running bands.x as:    mpirun -np 4 /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin/bands.x  -nk 1 -nd 1 -nb 1 -nt 1 
  running plotband.x as: /shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/bin/plotband.x

  running the scf calculation... done
  running pp.x to do a 2-d plot of the charge density... done
  running plotrho.x to generate rho.ps... done

  running pp.x to do another 2-d plot of the charge density... done
  generating si.charge.png...Could not find/open font when opening font "arial", using internal non-scalable font
 done
  generating contour plot of the charge si.contour.ps... done

  running the band-structure calculation for Si... done
  running the post-processing for band structure... done
  running plotband.x to generate sibands.ps... done
  creating the postscript file si.bands.ps... done

  cleaning /tmp/ucl-user_tempdir...
/home/ucl-dir/test/q-e/PP/examples/example01: done
```

**Congratulations! you have completed your first QE job! **

you can then enter the results folder and inspect the outcome. since there are some postscript and files, it may be best if you copy the entire folder to your pc, inspect the results there (using scp, etc, as shown in [file_howto.md](files_howto.md). 


# Running your own scripts
Inside the results folder from example 01 you will see that the script has generated some stand alone input files for QE. The first one we will look at is `si.scf.in` which is the simple self consistent field calculation for Silicon. Take a look at the file, and learn the various structure using the documentation. This will be the basis for your own custom files. We will do this in a seriest of tutorials. 
TIP: you can view a file by printing it and scrolling through it with the arrow keys using the `more` command like so: 

```sh 
(aristotle01 results) > more si.scf.in 
 &control
    calculation='scf'
    restart_mode='from_scratch',
    prefix='si'
    pseudo_dir = '/shared/ucl/apps/quantum-espresso-thermo_pw/6.5-fixedlinks/intel-2018/q-e/pseudo/',
    outdir='/tmp/ucl-user_tempdir/'
 /
 &system
    ibrav= 2, celldm(1)= 10.2, nat= 2, ntyp= 1,
    ecutwfc =18.0
 /
 &electrons
    conv_thr =  1.0d-8
    mixing_beta = 0.7
 /
ATOMIC_SPECIES
 Si  28.086  Si.pz-vbc.UPF
ATOMIC_POSITIONS alat
 Si 0.00 0.00 0.00
 Si 0.25 0.25 0.25
K_POINTS 
  10
   0.1250000  0.1250000  0.1250000   1.00
   0.1250000  0.1250000  0.3750000   3.00
   0.1250000  0.1250000  0.6250000   3.00
   0.1250000  0.1250000  0.8750000   3.00
   0.1250000  0.3750000  0.3750000   3.00
   0.1250000  0.3750000  0.6250000   6.00
   0.1250000  0.3750000  0.8750000   6.00
   0.1250000  0.6250000  0.6250000   3.00
   0.3750000  0.3750000  0.3750000   1.00
   0.3750000  0.3750000  0.6250000   3.00
(aristotle01 results) >
```

to run a stand alone QE input file do the following (see also [lammps getting started guide](lammps_getting_started.md) for how to run jobs in the foreground and background and the meaning of '>' and '<' shell operators)


```sh
pw.x < si.scf.in > si.scf.out
```

if you inspect the out file you can find for example the total energy of the configuration. 

## Exercise 
* If you change the position of the Si atoms slightly, would the energy increase or decrease? try it. 
* TIP: in order to get the total energy easilly without opening the file by hand, we can use the `grep` command like so:
```sh
(aristotle01 results) >     grep 'total energy' si.scf.out
     total energy              =     -15.84096534 Ry
     total energy              =     -15.84406358 Ry
     total energy              =     -15.84451017 Ry
     total energy              =     -15.84452626 Ry
     total energy              =     -15.84452724 Ry
!    total energy              =     -15.84452726 Ry
     The total energy is the sum of the following terms:
```
you see in the above that there are 6 self consistent field iterations, and the final one is the actual result for the energy. 
