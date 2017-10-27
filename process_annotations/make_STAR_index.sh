#!/bin/bash

#The job should run on the testing partition
#SBATCH -p main

#The name of the job is test_job
#SBATCH -J make_STAR_index

#The job requires 1 compute node
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=40000

#The maximum walltime of the job is a half hour
#SBATCH -t 07:59:00

#These commands are run on one of the nodes allocated to the job (batch node)
module load star-2.5.2
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir annotations/GRCh38/STAR_index_v90_oh74/ --genomeFastaFiles annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --sjdbGTFfile annotations/GRCh38/Ensembl_90/Homo_sapiens.GRCh38.90.gtf --sjdbOverhang 74

