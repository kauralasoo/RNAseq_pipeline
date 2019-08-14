#!/bin/bash

#The job should run on the testing partition
#SBATCH -p main

#The name of the job is test_job
#SBATCH -J run_susie

#The job requires 1 compute node
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 36:00:00

#SBATCH --mem 12000

#These commands are run on one of the nodes allocated to the job (batch node)
module load singularity
singularity exec susie-finemapping.img Rscript susie_finemapping.R

