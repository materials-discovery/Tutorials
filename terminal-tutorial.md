# First Stpes

- Open a terminal (e.g. using the materials-discovery.org) and 
type the following (following each command, press return or enter to execute)

```bash
# this will always take you back to the root, i.e., your home folder

$ cd 

# check which folder you are in

$ pwd
# this should be e.g., 
/home/jovyan

# list all files 
ls -larthF

# add this as an alias to .bashrc (remember this is all in the home folder!)
nano .bashrc

# look for the line with aliases, and add/change the existing ll
alias ll='ls -lrthaF'

# exit nano (pressing Ctr X, and then answering the questions on the screen)
# to make the changes active, run 
source .bashrc

```


- Moving, copying and deleting files
```bash
# again go to the home folder 
cd
touch eg.py
ll
less eg.py # this shows the contents. 
nano eg.py


```

- now lets create our first python script 
```bash
nano eg.py
# now add the following line:
```

```python 
print("Hello World!")
```

```bash
# now run in the terminal like so
python3 eg.py
# you should get Hello World! printed 
# next modify the program to type your name

```
```python 
myname="This is myname!"
print(f"Hello World!, my name is {myname}")
```


- lets create a folder and move our `eg.py` to it
```bash
# make sure you are in the root folder (or home) 
pwd
mkdir noname
ll # or ls 
# since this is a noname, lets remove it, 
rmdir noname 
ll
# lets try again 
mkdir alpha_fold
ll
# lets rename it (or really, move it) 
mv alpha_fold example_1
ll
# then cd into example_1
cd example_1
pwd

# lets move eg.py to example_1
mv eg.py example_1
ll
# now make a copy of eg.py
cp eg.py eg_copy.py
ll

# lets remove the copy

rm eg_copy.py
ll

```
## Important! change the default behavious of `rm, cp, and mv` to be more safe
- open the .bashrc (remember it is in your home folder) in nano or any other text or code editor. 
- search for teh alias section
- add the following or change the aliases as needed:
```bash
alias rm='rm -i'
alias cp='cp -i'
alias mv='mv -i'
```
## browsing folders 
- to go one folder up, type the command `cd ..`
- to go two folders up, `cd ../../`
- to go to the last visited folder do `cd -`

