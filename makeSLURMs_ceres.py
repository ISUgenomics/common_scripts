#!/usr/bin/python

Usage = """
creates a job file for desired number of commands per job

Usage:

  makeSLURMs.py <number of jobs per PBS file> <commands file>

eg:
  
  makeSLURMs.py 10 bowtie2.cmds

will create bowtie2_N.sub files, where N equals to number of lines in bowtie2.cmds divided by 10

If you have large number of commands that you would like to package (a set number) in a single
PBS script file, you can run this script along with desired number of commands per job.
Note that all commands will run in serially with this script (s suffix). If you want to run all commands at a time,
parallel fashion, then use the p suffix script.


PLEASE EDIT THIS SCRIPT TO INCLUDE THE RIGHT EMAIL! 

Arun Seetharam
arnstrm@iastate.edu
11/08/2016
"""
import sys
import os
if len(sys.argv)<3:
    print Usage
else:
   cmdargs = str(sys.argv)
   cmds = open(sys.argv[2],'r')
   jobname = str(os.path.splitext(sys.argv[2])[0])
   filecount = 0
   numcmds = int(sys.argv[1])
   line = cmds.readline()
   while line:
        cmd = []
        while len(cmd) != int(sys.argv[1]):
                cmd.append(line)
                line = cmds.readline()
        w = open(jobname+'_'+str(filecount)+'.sub','w')
        w.write("#!/bin/bash\n")
        w.write("#SBATCH --partition=short\n")
        w.write("#SBATCH --mem-per-cpu=3G\n")
        w.write("#SBATCH -n 20\n")
        w.write("#SBATCH -t 48:00:00\n")
        w.write("#SBATCH -J "+jobname+"_"+str(filecount)+"\n")
        w.write("#SBATCH -o "+jobname+"_"+str(filecount)+".o%j\n")
        w.write("#SBATCH -e "+jobname+"_"+str(filecount)+".e%j\n")
        w.write("#SBATCH --mail-user=email@gmail.com\n")
        w.write("#SBATCH --mail-type=BEGIN,END,FAIL\n")
        w.write("cd $SLURM_SUBMIT_DIR\n")
        w.write("ulimit -s unlimited\n")
        count = 0
        while (count < numcmds):
           w.write(cmd[count])
           count = count + 1
        w.write("scontrol show job $SLURM_JOB_ID\n")
        w.close()
        filecount += 1
