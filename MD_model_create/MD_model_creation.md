For molecular dynamics (MD) simulations we need to create a model system containing the atoms and  
- their types 
- and initial positions
- and velocities
- as well as define the simulation box 
- and the boundary conditions 

Usually you want to create the model in a format that your intended MD code package can understand. Note, in the near future we will switch to using the SimPhoNy package developed by our group and Fraunhofer, which uses a generalised standard model, but for now we focus on LAMMPS supported models. 

Take a look at the LAMMPS data format here: https://docs.lammps.org/99/data_format.html

In this tutorial we will do the following: 
1. Use the http://www.crystallography.net/ database to search and download some materials models in CIF format (Crystallographic Information File is the standard format for storing crystallographic structural data, see https://www.ccdc.cam.ac.uk/community/access-deposit-structures/deposit-a-structure/guide-to-cifs/)
2. we will specifically download a simple metal (Cu) CIF file (any one will do as long as the space group is Fm-3m)
	1. we are going to use this one for example: http://www.crystallography.net/cod/4105040.cif
3. If you are using your own computer, please download and install [VESTA](VESTA.md) as in this tutorial [VESTA Tutorial](VESTA.md). If you are using UCL@Anywhere desktop/machine, VESTA is already available, so no need to reinstall. 
4. Open VESTA, then load the Cu CIF file and visualise it
Below are some images of the steps: 

![[Pasted image 20230118112755.png]]
The VESTA window after loading the above Cu CIF file.
5. We will then duplicate the file and export it. Click on `Edit->Edit Data->Unit Cell-> Transform`  (note in the image below we changed the default of the Cu atoms, you can leave it as is though)
![[Pasted image 20230118114415.png]]
6. change the Rotation matrix so the diagonal elements are all equal to 4 (this means we will duplicate the cell by translating it along the original conventional unit cell latter directions 4 times).
	-  Quiz: check the description in the window for what the transformation is doing, do you know why it actually makes a translation in this case? 
![[Pasted image 20230118114701.png]]

7. click ok and apply, this is what you should get (note the dot dashed lines indicate the "unit cell" dimensions hace changed, this is what we call a [Super Cell](supercell.md))
![[Pasted image 20230118114731.png]]

1. we will save it and then use an application called ATOMSK to convert it to LAMMPS data file
	1. Go to `File -> Save as`, and save it to your [working folder](working folder) as a VESTA File, useful in case you want to open it later using VESTA and make additional changes to it, or you would need to recreate it from the start) 
	2. Go to `File -> Export Data`  and choose CIF file in the format box (see image below)  and save to the working folder. 
	![[Pasted image 20230118120356.png]]
1. let us convert the exported file to LAMMPS datafile using the tool [ATOMSK](https://atomsk.univ-lille.fr) which is open source  code developed by  Pierre Hirel.  
	1. now to use ATOMSK we need to move the file to [aristotle.rc.ucl.ac.uk]() or if you are using your laptop you can follow [ATOMSK Tutorial](ATOMSK.md) to install it and run it. If you need to work with ARISTOTLE follow this link [moving a file to ARISTOTLE] then 
	2. either go to ARISTOTLE or on you machine pen a terminal (in either windows or MAC or linux) and perform (type) the following:
	`% atomsk 4105040.cif 4105040.lmp `
	![[Pasted image 20230118122117.png]]
	Now you should have a file called `4105040.lmp` 
	6. use the `unix/linux` command `more` or `less` to view it: 
	![[Pasted image 20230118122240.png]]
	inspect  the file: [4105040.lmp](4105040.lmp) and compare to the [lammps manual](https://docs.lammps.org/99/data_format.html)
	8. by the way, ATOMSK can do much more, you can explore it here: [Explore ATOMSK]()

4. load it to VMD and Visualise it 
5. load it to OVITO and visualise it 
6. now we will move to [Optimise lattice parameter for CU](linux_tutorials/MD_tutorial_1/Cu.md) tutorial
	1. we will also look at some available tutorials that use instead the LAMMPS built in crystal maker functionality. [WIP](Tchupp_tutorial.md)






# Advanced topics
these are some typical and common topics, however for the scope of our moudle we refer to them as more advanced, and will be used mainly by those involved in MD projects. 
## Create a rotated unit cell for surface/interface models
this is useful for example if you are interested in free surfaces or interfaces with specific orientation relationships, e.g., a Cu (111) surface would require a different unit cell orientation that the standard conventional one you probably downloaded from the COD database above. 

## Create an Interface by joining different models 
 
