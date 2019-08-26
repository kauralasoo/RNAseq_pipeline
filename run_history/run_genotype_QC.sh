#ROSMAP illumina
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/illumina_raw/chop.rosmap.euam.vFinal.382only\
 --output_name ROSMAP_illumina_GRCh37_genotyped\
 --outdir ROSMAP_illumina

#ROSMAP affy
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/plink_raw/ROSMAP_raw_GRCh37_bed_filtered\
 --output_name ROSMAP_affy_GRCh37_genotyped\
 --outdir ROSMAP_affy

# Schmiedel_2018
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/Schmiedel_2018/genotypes/matrix/IC_DNA\
 --output_name Schmiedel_2018_GRCh37_genotyped\
 --outdir Schmiedel_2018

# Fairfax_2018
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/Fairfax_2018/genotyes/raw_genotyped/Fairfax_2018\
 --output_name Fairfax_2018_GRCh37_genotyped\
 --outdir Fairfax_2018

 # BrainSeq h650
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/BrainSeq/genotypes/processed/BrainSeq_h650\
 --output_name BrainSeq_h650_GRCh37_genotyped\
 --outdir BrainSeq_h650 

# BrainSeq 1M
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/BrainSeq/genotypes/processed/BrainSeq_1M\
 --output_name BrainSeq_1M_GRCh37_genotyped\
 --outdir BrainSeq_1M 


#### Post-imputation
# ROSMAP_illumina
nextflow run crossmap_genotypes.nf -profile crossmap -resume\
 --vcf_files "/gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/Michigan_GRCh37_Phase3_200819/illumina/GRCh37/chr*.dose.vcf.gz"\
 --output_name ROSMAP_illumina 

# BrainSeq 1M
nextflow run crossmap_genotypes.nf -profile crossmap -resume\
 --vcf_files "/gpfs/hpc/home/a72094/datasets/controlled_access/BrainSeq/genotypes/Michigan_GRCh37_Phase3_250819/1M/GRCh37/chr*.dose.vcf.gz"\
 --output_name BrainSeq_1M\
 --outdir BrainSeq_1M

# ROSMAP_affy
 nextflow run crossmap_genotypes.nf -profile crossmap -resume\
 --vcf_files "/gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/Michigan_GRCh37_Phase3_200819/affy/GRCh37/chr*.dose.vcf.gz"\
 --output_name ROSMAP_affy\
 --outdir ROSMAP_affy
