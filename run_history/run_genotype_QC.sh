#ROSMAP illumina
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/illumina_raw/chop.rosmap.euam.vFinal.382only\
 --output_name ROSMAP_illumina_GRCh37_genotyped\
 --outdir ROSMAP_illumina

#ROSMAP affy
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/ROSMAP/genotypes/plink_raw/ROSMAP_raw_GRCh37_bed\
 --output_name ROSMAP_affy_GRCh37_genotyped\
 --outdir ROSMAP_affy

# Schmiedel_2018
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/Schmiedel_2018/genotypes/matrix/IC_DNA\
 --output_name Schmiedel_2018_GRCh37_genotyped\
 --outdir Schmiedel_2018
