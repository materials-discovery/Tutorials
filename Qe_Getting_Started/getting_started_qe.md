# Quantum ESPRESSO on Materials Discovery Platform

# Introduction 
This is a getting started tutorial/howto for setting up and running quantum Espresso [https://www.quantum-espresso.org](https://www.quantum-espresso.org) on Materials Discovery. 

## What is QE
* Quantum ESPRESSO official website:  [https://www.quantum-espresso.org](https://www.quantum-espresso.org). 

![https://www.quantum-espresso.org](/figures/quantum_ogo_ok.png)

> Quantum ESPRESSO is an integrated suite of Open-Source computer codes for electronic-structure calculations and materials modeling at the nanoscale. It is based on density-functional theory, plane waves, and pseudopotentials.



## Where can I find more information about QE
The best source of information is the [Quantum ESPRESSO home page](https://www.quantum-espresso.org/) and the [documentation](https://www.quantum-espresso.org/documentation/). Quantum ESPRESSO is also found on the official [github repository](https://github.com/QEF/q-e).

Also see the official documentation and tutorials page: [https://www.quantum-espresso.org/tutorials/](https://www.quantum-espresso.org/tutorials/)

## Prerequisite 
In order to follow the tutorials here you should first read the `terminal-tutorial.md` in this repo and get aquatinted with using the terminal and linux shell. 

# Check Quantum Espresso available on Materials discovery
Type in the command line the following for the plane wave pseudo potential code:

```shell
which pw.x 
```
If you see something like this, you are all set for simple qe jobs:

```shell
(base) jovyan@c23b95e7ecee:~$ which pw.x
/usr/bin/pw.x
(base) jovyan@c23b95e7ecee:~$ 
```

# Run the first example

## Download Standard solid-state pseudopotentials (SSSP)
Visit [Materials Cloud SSSP Efficiency page](https://www.materialscloud.org/discover/sssp/table/efficiency) to download the entire pseudopotential folder (e.g., SSSP_1.3.0_PBE_efficiency). After downloading, upload the required files to the Materials Discovery platform for further use. Ensure that the pseudopotentials match the calculation settings accurately. In this case, we upload pseudopotential file of `Al`,`Cu`,`Si`. Please copy the path of your pseudopotential file using `pwd` command. 


### View and modify the input file
You can use `nano` command like below or `vim`
```
nano si.scf.in
```
Modify the pseudopotential path and filename in the input file to match the location of the uploaded files. For example, locate the `PSEUDO_DIR variable` or the `pseudopotential file reference` (e.g., Al.pbe-n-kjpaw_psl.1.0.0.UPF) and update it with the correct path. If the pseudopotentials are stored in the folder SSSP_1.3.0_PBE_efficiency, you might update it like this: `PSEUDO_DIR ="/path/to/SSSP_1.3.0_PBE_efficiency"`, Please paste the path with your own.

Also, modify the pseudopotential filename in the input file `si.scf.in`, using `nano` or `vim`. if your filename is `Si.pbe-n-rrkjus_psl.1.0.0.UPF`, you should modify the file like below
```
ATOMIC_SPECIES
 Si  28.086  Si.pbe-n-rrkjus_psl.1.0.0.UPF
```

Save it.

## Run the input file
Make sure you are in `Qe_Getting_Started` folder

Run a stand alone QE input file do the following (see also lammps getting started guide for how to run jobs in the foreground and background and the meaning of '>' and '<' shell operators)
```
pw.x < si.scf.in > si.scf.out
```
# Congrutulations, you finish your first DFT calculation. 

# Dig into the input, and output file
![input file](/figures/input.jpeg)
![output file](/figures/output.png)
# Run Your own Script
You can find more examples on Quantum ESPRESSO  official [github repository](https://github.com/QEF/q-e/tree/develop/PW/examples). And please make sure you have the right pseudopotential path and filename in the input file. Also, check the parameters you don't understand on [Quantum espresso website](https://www.quantum-espresso.org/Doc/INPUT_PW.html)
