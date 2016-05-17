#!/usr/bin/env zsh

### Job name
#BSUB -J LEVELSET

### File / path where STDOUT will be written , the %J is the job id
#BSUB -o LEVELSET.%J

### (OFF ) Different file for STDERR , if not to be merged with STDOUT
#BSUB -e LEVELSET.e%J



###  Request vitual memory  you  need  for your job  in  MB / 50000 = 9GB

#BSUB -M 9000


### Request the time you need for execution in minutes
### The format for the parameter is: [ hour :] minute ,
### that means for 80 minutes you could also use this : 1:20
#BSUB -W 3:00

### Specify your mail address
#BSUB -u miessen@imm.rwth-aachen.de

### Send a mail when job is done
#BSUB -N

#### compute nodes
#BSUB -n 8


### for shared memory jobs (OpenMP)
#BSUB -a openmp

### Change to the work directory
mkdir mySim_$J 
cd mySim_$J
cp ../LevelSet_IMM ./
cp ../parameters.xml ./

### Execute your application
numactl --interleave=all LevelSet_IMM parameters.xml
