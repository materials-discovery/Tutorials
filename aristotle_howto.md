# Aristotle@UCL

This is a tutorial for working with the Aristotle computer service (read: a server, a remote computer) on the UCL computer center resources. 

## What is Aristotle
Aristotle is an interactive, Linux-based compute service (a computer) for teaching, running on three each node 64 gigabytes of RAM and 16 cores. It is maintained by the UCL compute services and is accessible for all students and teaching staff. 

The nodes run the Red Hat Enterprise Linux operating system (RHEL 7). 

For more details see the Information Services Division web pages at UCL: [ISD](
 https://www.rc.ucl.ac.uk/docs/Other_Services/Aristotle/) for more information.

## How to connect to Aristotle from your PC/MAC/LINUX machine

for the most up to date information, see the official documentation on [UCL-IDS website](https://www.rc.ucl.ac.uk/docs/howto/). An excerpt is below, but always refer or at least check the official docs. 

## If you are using a Linux / Unix / Mac OS X / Windows 10 or 11

Use the terminal application and type the below command to secure shell (ssh) into the machine you wish to access. 

Replace <your_UCL_user_id> with your central UCL username, and <system_name> with the name of the machine you want to log in to, eg. myriad, kathleen, Aristotle.


`ssh <your_UCL_user_id>@<system_name>.rc.ucl.ac.uk`

**NOTE**: if you are not connected to the EDUROAM wifi network while at UCL, then you need to use [VPN](https://www.ucl.ac.uk/isd/services/get-connected/ucl-virtual-private-network-vpn) first before connecting to ARISTOTLE (alternatively you van use a gateway but is outside the scope of this how-to)

**NOTE** for connecting from China, please see the official [UCL documentation](https://www.rc.ucl.ac.uk/docs/howto/#china-connect)

**You need to use a Terminal for running applications on Aristotle!**

For this tutorial you need a Terminal program running an appropriate shell on your laptop/desktop and connect to a shell (a remote terminal) in Aristotle. ![](tutorials/Tutorials/figures/Linux_Terminal.png).

This is also needed if you are using [UCL Desktop@Anywhere](https://www.ucl.ac.uk/isd/services/computers/remote-access/desktopucl-anywhere)

### Additional information for Microsoft Windows PC users

On Windows you need something that will give you a suitable terminal and ssh - usually PuTTY, or on Windows 10/11 you can use the builtin OpenSSH from a command prompt and type the same ssh command as the Linux instructions.

* Get [putty](https://www.rc.ucl.ac.uk/docs/howto/#windows-putty). For details how to set it up see: [https://www.rc.ucl.ac.uk/docs/howto/](https://www.rc.ucl.ac.uk/docs/howto/)
* or on windows 100/11 open a windows "CMD" or Power-shell. To use CMD, go to the search area in the panel and type cmd or command, it should appear like so: 

![](tutorials/Tutorials/figures/cmd_search2.png) 

then click on Command Prompt to launch it, and type in it 

`ssh <ucl_user_id>@aristotle.rc.ucl.ac.uk`

(you can also use powershell instead). 

### Additional information for MacOs and Linux users

For **Mac/Linux** open a Terminal application (use spotlight search to find it) then use ssh to connect to Aristotle as in windows CMD above, i.e., type:

`ssh <ucl_user_id>@aristotle.rc.ucl.ac.uk`

### Enter your UCL password once asked

Then you will be asked for the password. **Note while you type the password on the terminal, you will not see anything being typed! this is to make sure now one sees your password! just type the password carefully and then hit Enter key**

Now you should be logged in and get something similar to the following prompt: 

![](tutorials/Tutorials/figures/cmd.png)

### More information about using a Unix(R) like shell - The Command Line Interface (CLI)
When you log into Aristotle you are actually getting a terminal with a command line interface (CLI). A CLI is a text-based interface that allows users to manage the computer and run applications, and organize files. A CLI expects an input command, executes it, and prints the output. Specifically, once you ssh into Aristotle you will get a command line with a shell CLI  running in it, usually this will be a bash shell. 

For much more thorough information and help about the linux/Mac (in fact, more generally UNIX) Terminal see [A guide to the Linux terminal for beginners](https://opensource.com/article/21/8/linux-terminal) and a nice intro to the [main commands you need in a terminal](https://www.redhat.com/sysadmin/10-commands-terminal).

# Copying files from/to your PC/MAC to/from Aristotle 

## Using a terminal
similar to the `ssh` package (secure shell) for logging in, there is a secure copy package called `scp` which allows one to copy files between machines using a terminal. It is invoked similar to ssh. 

In general (commands are always run in a Terminal when preceded with `>$ ` )

`scp -r <PATH-TO-WINDOWS-FILE> <ucl-user>@aristotle.rc.ucl.ac.uk:/path/to/where/you/want/the/file`

where `<PATH-TO-WINDOWS-FILE> ` should be replaced with the current windows path (can be copied from the path in file-explorer) and `<ucl-user>` should be replaced by your ucl user name. 

note on unix the folder separator is a forward slash (`/`) unlike in windows which is a backslash (`\`) 

For example suppose you want to copy `somefile.in` in your home folder in windows to a folder called ProjectX on Aristotle:

`scp -r C:\Users\username\somefile.in  <ucl-user>@aristotle.rc.ucl.ac.uk:/home/<ucl-user>/ProjectX/`

Note that `/home/<ucl-user>` is your home user folder on Aristotle (and in, practically, all unix systems)
to copy the file back from Aristotle to your PC, just reverse, like so: 

`scp -r C:\Users\username\somefile.in  <ucl-user>@aristotle.rc.ucl.ac.uk:/home/<ucl-user>/ProjectX/`


## Using a GUI programme

you can use a few GUI (graphical user interface) programs like winSCP [https://winscp.net/eng/index.php](https://winscp.net/eng/index.php)
or MobaXterm [https://mobaxterm.mobatek.net](https://mobaxterm.mobatek.net) on windows or Filezilla [https://filezilla-project.org/](https://filezilla-project.org/) which works also on a mac and linux. See the respective pages for how to, and ask the PGT in case you need further help.

## Other information 
see [https://www.ucl.ac.uk/isd/accessing-research-data-storage-through-ssh-windows](https://www.ucl.ac.uk/isd/accessing-research-data-storage-through-ssh-windows)

## Setting up and running Lammps on Aristotle 
See the file `aristotle_lammps_howto.md` in this repo. 

## Setting up and running Quantum Espresso on Aristotle 
See the file `aristotle_quantumEspresso_howto.md` in this repo. 


# Loading standard modules into the Aristotle environment 

Since Aristotle has to accommodate a large number of packages and programs each having often different dependencies on other packages or software libraries, Aristotle utilizes a program called `module` to manage the applications available to the user. 
A standardized set of packages and of useful tools should be loaded first including Intel compilers but users can also change their own environment.

In order to make life easier, load the following initial default modules by the runing the command 

``` module load default-modules ```

inside your terminal (the one you ssh-ed into above) 

You could also add the above line into you `.bashrc` file which is found in your home folder (see below how to do this) so it is loaded on startup (actually at login).

with this command loaded, running the command (in the terminal)

```module list ```

should give something like  

```
Currently Loaded Modulefiles:
  1) ops-tools/2.0.0               12) nedit/5.6-aug15
  2) gcc-libs/4.9.2                13) dos2unix/7.3
  3) cmake/3.21.1                  14) giflib/5.1.1
  4) flex/2.5.39                   15) emacs/26.3
  5) git/2.32.0                    16) tmux/3.2a
  6) apr/1.7.0                     17) mrxvt/0.5.4
  7) apr-util/1.6.1                18) userscripts/1.4.0
  8) subversion/1.14.1             19) rcps-core/1.0.0
  9) screen/4.8.0-ucl1             20) compilers/intel/2018/update3
 10) gerun                         21) mpi/intel/2018/update3/intel
 11) nano/2.4.2                    22) default-modules/2018
```

# Other useful information 

a quick search on youtube can return many useful demonstrations of using the terminal and the bash shell. 

more coming soon.. 