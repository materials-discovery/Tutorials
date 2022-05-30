# General commands for Git

* Note: Terminal refers to a Linux/Unix/ shell terminal. 
* All commands below are supposed to be run in a Unix (bash or zsh) shell in your machine or remotely. 

|command  | description   |
|--|--|
| `git remote -v -a`                    | Run in a Terminal in local git repository. Shows the connected remote (usually on github) |
| `alias gitdiff="git difftool --tool=vimdiff ` | add this to your .aliases or .bashrc (in aristotle)/.zshrc (on mac) to see a diff side by side. You can also run the command directly in a terminal so: `  git difftool --tool=vimdiff file1 file2` etc. |
| `git add -u ` | update all changes done in local directory instead of doing `git add file1 file2...`, useful especially when one moves or deletes a number of files|
| `git add newFile` | if you create a new source code file, then you need to add the new file it to thr git system (done locally in your own repository) |
| `git fetch --all` | git fetch all branches from the remote (GitHub usually in our case) so that subsequent git checkout remote-branch-name will actually get the remote branch and not create a new one, so always do this before checking (switching) to a remote branch. |
| `git checkout branch-name`| switch to the existing branch (always do git fetch --all before hand) from remote (GitHub). usually should be followed by git pull |
| `git pull` | pull the current branch from the remote |
 

## coming soon: 
- how to install git-bash under windows
- how to install and use WSL2 under windows. 
- how to use jupyter notebook...  