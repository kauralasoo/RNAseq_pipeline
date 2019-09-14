 ##### Naranbhai_2015 ####
nextflow run main.nf -profile eqtl_catalogue\
 --expression_matrix /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/Naranbhai_2015.HumanHT-12_V4_norm_exprs.tsv\
 --phenotype_metadata /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 --sample_metadata /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv\
 --genotype_vcf /gpfs/hpchome/a72094/datasets/controlled_access/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
 --is_imputed true\
 -resume

  ##### Fairfax_2014 ####
nextflow run main.nf -profile eqtl_catalogue\
 --expression_matrix /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/Fairfax_2014.HumanHT-12_V4_norm_exprs.tsv\
 --phenotype_metadata /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 --sample_metadata /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Fairfax_2014.tsv\
 --genotype_vcf /gpfs/hpchome/a72094/datasets/controlled_access/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
 --is_imputed true\
 -resume

##### Fairfax_2012 ####
nextflow run main.nf -profile eqtl_catalogue\
 --expression_matrix /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/Fairfax_2012.HumanHT-12_V4_norm_exprs.tsv\
 --phenotype_metadata /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 --sample_metadata /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Fairfax_2012.tsv\
 --genotype_vcf /gpfs/hpchome/a72094/datasets/controlled_access/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
 --is_imputed true\
 -resume

##### CEDAR ####
nextflow run main.nf -profile eqtl_catalogue\
 --expression_matrix /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/CEDAR.HumanHT-12_V4_norm_exprs.tsv\
 --phenotype_metadata /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 --sample_metadata /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/CEDAR.tsv\
 --genotype_vcf /gpfs/hpchome/a72094/datasets/open_access/CEDAR/genotypes/Michigan_GRCh37_1KGPhase3_220918/GRCh38/CEDAR_GRCh38.filtered.renamed.vcf.gz\
 --is_imputed true\
 -resume

##### Kasela_2017 ####
nextflow run main.nf -profile eqtl_catalogue\
 --expression_matrix /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/Kasela_2017.HumanHT-12_V4_norm_exprs.tsv\
 --phenotype_metadata /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 --sample_metadata /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Kasela_2017.tsv\
 --genotype_vcf /gpfs/hpc/home/a72094/datasets/controlled_access/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.filtered.vcf.gz\
 --is_imputed true\
 -resume


 ##### GENCORD ####
 ## transcript_usage
nextflow run main.nf -profile eqtl_catalogue\
 --expression_matrix /gpfs/hpchome/a72094/datasets/processed/transcript_usage/GENCORD.transcript_usage_qnorm.tsv\
 --phenotype_metadata /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz\
 --sample_metadata /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 --genotype_vcf /gpfs/hpc/home/a72094/datasets/controlled_access/GENCORD2/genotypes/Michigan_GRCh37_1KGPhase3_220918/GRCh38/GENCORD_GRCh38.filtered.vcf.gz\
 --is_imputed true\
 -resume

 #Run multi_study version
 #HipSci
nextflow run main_multi_study.nf -profile eqtl_catalogue\
 --readPathsFile /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/qtlmap/multi_study_run_Hipsci.tsv\
  --outdir summary_statistics\
  --n_batches 400\
  --run_permutation false\
  -resume

#Ye_2018
nextflow run main_multi_study.nf -profile eqtl_catalogue\
 --readPathsFile /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/qtlmap/multi_study_run_Ye_2018.tsv\
 --outdir summary_statistics\
 --n_batches 400\
 --run_permutation false\
 -resume


nextflow run main_multi_study.nf -profile eqtl_catalogue --readPathsFile /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/qtlmap/multi_study_run_Hipsci.tsv --outdir summary_statistics --n_batches 400 --run_permutation true --run_nominal false -resume

#Quach and Schmiedel permutation runs
nextflow run main_multi_study.nf -profile eqtl_catalogue\
 --studyFile /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/qtlmap/multi_study_run_big.tsv\
 --outdir summary_statistics\
 --n_batches 400\
 --is_imputed true\
 --run_permutation true\
 --run_nominal false\
 -resume