See: https://atomsk.univ-lille.fr/
*Atomsk is a free, Open Source command-line program dedicated to the creation, manipulation, and conversion of data files for atomic-scale simulations in the field of computational materials sciences. It supports the file formats used by many programs, including:

*Visualization and analysis softwares: Atomeye, ddplot, OVITO, VESTA, XCrySDen;
Ab initio calculation softwares: ABINIT, Quantum Espresso, VASP, SIESTA;
Classical force-field simulation softwares: DL_POLY, GULP, IMD, LAMMPS, XMD;
TEM image simulation softwares: Dr Probe, JEMS, QSTEM.

![[Pasted image 20230118185648.png]]


- First lets install it on Aristotle 
1. open a terminal to Aristotle and connect in addition with winscp (if on mac you can use )
2. go to the github repository and clone the repository locally like so: 
```shell
git clone git@github.com:pierrehirel/atomsk.git
```
![[Pasted image 20230118190144.png]]

3. now cd into the new folder `% cd atomsk`
4. you can view the `README` file using winscp viewer or the terminal using `less README ` (remember to click q to exit `less` when you finish)
5. according to the readme file, to install we need to enter the src folder and do `make atomsk` this will run the gfortran compiler and compile the code producing the executable atomsk which we will add to the path like so: 
6. add the line `export PATH=$PATH:$HOME/atomsk/src` (replace ucsaha with your user name! ) to the end of the file ` .bashrc` which is in your home folder. (hint: `nanao ~/.bashrc`) then save and run the following in the terminal `% . ~/.bashrc` to activate the change. now if you type `%atomsk -h` you will get a message referring you to the atomsk website for information, this is a sign that the code is running well. 
![[Pasted image 20230118190719.png]]
