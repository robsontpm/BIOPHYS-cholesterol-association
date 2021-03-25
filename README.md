A set of scripts for computation of various parametyers in molecular dynamics 
simulation of the spontanous association of cholesterol and its products in water. 


INSTALATION
===========

1. This progam uses older python 2.7. We do not plan to port it to python 3.+ as of now. 
   The code depends on the numpy library (but very limited number of base functions),
   so any version of the library should do. 

2. You need to install library DMG-alpha first. You can download it here: dmgalpha.scirsc.org. 
   This page contains instalation instructions for the library. Script 'dipole_angles.py' 
   uses gromacs (http://www.gromacs.org/) to generate some data from the input trajectories. 
   Without gromacs, it will not work (only this script).


RUN CONFIGURATION
=================

3. The scripts use importlib to load setup for the run. THIS IS NOT A SAFE SOLUTION in open 
   environment, so for example you should not use it to create some kind of computational
   service availabe from web browser or similar. 

   This was deliberatelly used, as to have a simple way to write 'dynamic' 
   input files, where one could use arrays, generators, functions, etc. to
   write input files faster and not need to generate them with external
   scripts. 

4. Exemplary setup files are in 'sample_input'. You can prepare one file and
   it should be good for all scripts. There is a lot of documentation in the sample
   file and if you look at sample .pdb files you can get grasp what is going on. 
   It is advisable to use absolute paths, so that the script does not get confused. 
   The output directories will be created with a command equivalent to 'mkdir -p /path/to/dir'
   so they do not need to exist before running the scripts. By default the scripts always
   make a copy of the output folder if it exists. You can supply 'override' option
   to not make a copy. 

   Standard way to call a script:

   		python script_name.py path.to.conf.dir.filename

   (with backup) or 

   		python script_name.py path.to.conf.dir.filename override

   (without backup). The setup file in this example is called 'filename.py'
   and it is in the directory './path/to/conf/dir' (or any path that python
   recognizes as containing modules, in particular, there should be 
   '__init__.py' (can be empty) files in './path/to/conf/dir')
   

FOOTNOTES
=========

6. The scripts were written in linux and with intend to use it on that platform. 
   At the moment we do not plan to make a portable version for Mac or Windows. 




