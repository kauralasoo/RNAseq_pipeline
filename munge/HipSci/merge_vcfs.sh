#!/bin/bash

#The job should run on the testing partition
#SBATCH -p main

#The name of the job is test_job
#SBATCH -J bwa_index

#The job requires 1 compute node
#SBATCH -N 1

#The job requires 1 task per node
#SBATCH --ntasks-per-node=1

#The maximum walltime of the job is a half hour
#SBATCH -t 24:00:00

#SBATCH --mem 8000

#These commands are run on one of the nodes allocated to the job (batch node)
module load bcftools-1.9
bcftools merge open/chr1.dose.vcf.gz v1-0/chr1.dose.vcf.gz v1-1/chr1.dose.vcf.gz -Oz -o merged/chr1.dose.vcf.gz
bcftools merge open/chr2.dose.vcf.gz v1-0/chr2.dose.vcf.gz v1-1/chr2.dose.vcf.gz -Oz -o merged/chr2.dose.vcf.gz
bcftools merge open/chr3.dose.vcf.gz v1-0/chr3.dose.vcf.gz v1-1/chr3.dose.vcf.gz -Oz -o merged/chr3.dose.vcf.gz
bcftools merge open/chr4.dose.vcf.gz v1-0/chr4.dose.vcf.gz v1-1/chr4.dose.vcf.gz -Oz -o merged/chr4.dose.vcf.gz
bcftools merge open/chr5.dose.vcf.gz v1-0/chr5.dose.vcf.gz v1-1/chr5.dose.vcf.gz -Oz -o merged/chr5.dose.vcf.gz
bcftools merge open/chr6.dose.vcf.gz v1-0/chr6.dose.vcf.gz v1-1/chr6.dose.vcf.gz -Oz -o merged/chr6.dose.vcf.gz
bcftools merge open/chr7.dose.vcf.gz v1-0/chr7.dose.vcf.gz v1-1/chr7.dose.vcf.gz -Oz -o merged/chr7.dose.vcf.gz
bcftools merge open/chr8.dose.vcf.gz v1-0/chr8.dose.vcf.gz v1-1/chr8.dose.vcf.gz -Oz -o merged/chr8.dose.vcf.gz
bcftools merge open/chr9.dose.vcf.gz v1-0/chr9.dose.vcf.gz v1-1/chr9.dose.vcf.gz -Oz -o merged/chr9.dose.vcf.gz
bcftools merge open/chr10.dose.vcf.gz v1-0/chr10.dose.vcf.gz v1-1/chr10.dose.vcf.gz -Oz -o merged/chr10.dose.vcf.gz
bcftools merge open/chr11.dose.vcf.gz v1-0/chr11.dose.vcf.gz v1-1/chr11.dose.vcf.gz -Oz -o merged/chr11.dose.vcf.gz
bcftools merge open/chr12.dose.vcf.gz v1-0/chr12.dose.vcf.gz v1-1/chr12.dose.vcf.gz -Oz -o merged/chr12.dose.vcf.gz
bcftools merge open/chr13.dose.vcf.gz v1-0/chr13.dose.vcf.gz v1-1/chr13.dose.vcf.gz -Oz -o merged/chr13.dose.vcf.gz
bcftools merge open/chr14.dose.vcf.gz v1-0/chr14.dose.vcf.gz v1-1/chr14.dose.vcf.gz -Oz -o merged/chr14.dose.vcf.gz
bcftools merge open/chr15.dose.vcf.gz v1-0/chr15.dose.vcf.gz v1-1/chr15.dose.vcf.gz -Oz -o merged/chr15.dose.vcf.gz
bcftools merge open/chr16.dose.vcf.gz v1-0/chr16.dose.vcf.gz v1-1/chr16.dose.vcf.gz -Oz -o merged/chr16.dose.vcf.gz
bcftools merge open/chr17.dose.vcf.gz v1-0/chr17.dose.vcf.gz v1-1/chr17.dose.vcf.gz -Oz -o merged/chr17.dose.vcf.gz
bcftools merge open/chr18.dose.vcf.gz v1-0/chr18.dose.vcf.gz v1-1/chr18.dose.vcf.gz -Oz -o merged/chr18.dose.vcf.gz
bcftools merge open/chr19.dose.vcf.gz v1-0/chr19.dose.vcf.gz v1-1/chr19.dose.vcf.gz -Oz -o merged/chr19.dose.vcf.gz
bcftools merge open/chr20.dose.vcf.gz v1-0/chr20.dose.vcf.gz v1-1/chr20.dose.vcf.gz -Oz -o merged/chr20.dose.vcf.gz
bcftools merge open/chr21.dose.vcf.gz v1-0/chr21.dose.vcf.gz v1-1/chr21.dose.vcf.gz -Oz -o merged/chr21.dose.vcf.gz
bcftools merge open/chr22.dose.vcf.gz v1-0/chr22.dose.vcf.gz v1-1/chr22.dose.vcf.gz -Oz -o merged/chr22.dose.vcf.gz

