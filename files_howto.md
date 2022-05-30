### Inspect LAMMPS or Quantum ESPRESSO or other Input/output files on Aristotle

Say you need to inspect (or read, edit) `some-file` on Aristotle. Since this file is already on Aristotle we need some tools either on the machine itself, or we need to copy it to our desktop/laptop and look at it. 

One option is to use [scp or winscp]((https://www.rc.ucl.ac.uk/docs/howto/#windows-winscp)) and copy the file to your desktop/laptop, view, edit, change etc., using your windows tool which you are familiar with (you can use any tool which edits text, not wordpad or MS word as they will alter the format and make it void! then you can move it back to Aristotle. 

To edit such files on Windows, it is recommended to use a tool like [vscode](https://code.visualstudio.com) which is an open source solution from Microsoft, and can even do remote editing in addition to local files, see [https://code.visualstudio.com/docs/remote/ssh](https://code.visualstudio.com/docs/remote/ssh) for how to set it up for remote editing.

A simple solution is to directly view and edit the files on Aristotle i.e., without copying it over. For this we will use the simple text editor program called [nano](https://www.nano-editor.org) (alternatively, if you are familiar with `vim`, or are feeling adventures ;-), you may try to learn [vim](https://opensource.com/article/19/3/getting-started-vim) instead, which is already active by default). `nano` has to be activated first: 

* `module load nano`

now run 

* `nano some-file`  and use the arrow keys to move around.  
  

  Another Tip: instead of scp/winscp and a CMD terminal, you can use the more advanced, but easy to use  [MobaXterm](https://www.rc.ucl.ac.uk/docs/howto/#mobaxterm)) program!

