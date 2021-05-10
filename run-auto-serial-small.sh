#!/bin/bash

#First load all related modules.
#You can put the below two lines in a batch file.
#But remember the modules might get unloaded so check you loaded modules frequently.
module load intel/2018.4
module load intelmpi/2018.4

#Compile the code
make clean
make

#Schedule your jobs with sbatch
sbatch auto-teach-serial-small