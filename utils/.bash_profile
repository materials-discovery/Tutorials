# sumple .bash_profile, 
# Professor Adham Hashibon, a.hashibon@ucl.ac.uk

# If not running interactively, don't do anything
[ -z "$PS1" ] && return

  PS1="(\A \u@\h \W) > "

source ./.bashrc
