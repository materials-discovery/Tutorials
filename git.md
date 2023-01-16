## General commands for Git

* Note: Terminal refers to a Linux/Unix/ shell terminal. 
* All commands below are supposed to be run in a Unix (bash or zsh) shell in your machine or remotely. 

- `git clone url ` 
Clone a repository hosted at URL, e.g., 
```
git@github.com:materials-discovery/linux_tutorials.git 
```
This creates a copy of the folder on GitHub locally. 

- `git fetch ` 
Updates the metadata of the remote repo locally, so `git pull` gets the latest changes 
- `git fetch --all` 
git fetch all branches from the remote (GitHub usually in our case) so that subsequent `git checkout remote-branch-name` will actually get the remote branch and not create a new one, so always do this before checking (switching) to a remote branch. 
- `git pull` 
pull the current branch from the remote 
- `git remote -v `        
Run in a Terminal in local git repository folder. Shows the connected remote (usually on github) 
- `git add newFile` 
if you create a new source code file, then you need to add the new file it to thr git system (done locally in your own repository) 
- `git add -u ` | update all changes done in local directory instead of doing `git add file1 file2...`, useful especially when one moves or deletes a number of files
-  `alias gitdiff="git difftool --tool=vimdiff `
add this to your .aliases or .bashrc (in aristotle)/.zshrc (on mac) to see a diff side by side. You can also run the command directly in a terminal so: `  git difftool --tool=vimdiff file1 file2` etc. |
- `git checkout branch-name`
switch to the existing branch (always do git fetch --all before hand) from remote (GitHub). usually should be followed by git pull 


## Change origin
* Git push existing repo to a new and different remote repo server?

Create a new repo at github. Clone the repo from  to your local machine. Now you can work with it just like any other github repo. To pull in patches from upstream, simply run `git pull upstream master && git push origin master` [stackoverflow](https://stackoverflow.com/questions/5181845/git-push-existing-repo-to-a-new-and-different-remote-repo-server) |

```
git remote rename origin upstream 
git remote add origin URL_TO_GITHUB_REP
git push origin master
``` 
