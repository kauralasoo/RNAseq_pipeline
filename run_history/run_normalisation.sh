##### GENCORD ####
singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/GENCORD/featureCounts/merged_gene_counts.txt\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE

 ##### Lepik_2017 ####
 #gene_counts
singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/Lepik_2017/featureCounts/merged_gene_counts.txt\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Lepik_2017.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/gene_counts/\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE

 #transcript_usage
singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/Lepik_2017/Salmon/merged_counts/gencode.v30.transcripts/gencode.v30.transcripts.TPM.merged.txt\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Lepik_2017.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/transcript_usage/\
 -q transcript_usage\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE


 ##### Naranbhai_2015 ####
 singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/raw/Naranbhai_2015.tsv.gz\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/\
 -q HumanHT-12_V4\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE

##### Fairfax_2014 ####
 singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/raw/Fairfax_2014.tsv.gz\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Fairfax_2014.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/\
 -q HumanHT-12_V4\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE

##### Fairfax_2012 ####
 singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/raw/Fairfax_2012.tsv.gz\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Fairfax_2012.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/\
 -q HumanHT-12_V4\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE

##### Kasela_2017 ####
 singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/raw/Kasela_2017.tsv.gz\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Kasela_2017.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/\
 -q HumanHT-12_V4\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE

 ##### CEDAR ####
 singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/raw/CEDAR.tsv.gz\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/CEDAR.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/HumanHT-12_V4/\
 -q HumanHT-12_V4\
 --eqtlutils ../eQTLUtils/\
 --filter_qc TRUE
