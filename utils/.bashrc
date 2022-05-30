# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi


# User specific aliases and functions

force_color_prompt=yes 

export LS_COLORS=$LS_COLORS:'di=1;36:'
alias rm="/bin/rm -vi"
alias mv="/bin/mv -vi"
alias cp="/bin/cp -vi"
alias ls="ls -G --color"
alias ll="ls -salh"
alias l="ls -prth"



uncomment to load lammps in every log in
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


module load compilers/intel/2018/update3
module load mpi/intel/2018/update3
module load xorg-utils/X11R7.7
module load gnuplot/5.0.1
module load quantum-espresso/6.5-impi/thermo_pw-1.2.1/intel-2018


