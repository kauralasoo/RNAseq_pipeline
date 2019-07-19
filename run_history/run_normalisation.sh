##### GENCORD ####
singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/GENCORD/featureCounts/merged_gene_counts.txt\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/


 ##### Lepik_2017 ####
singularity exec qtl_norm_qc.img Rscript normaliseCountMatrix.R\
 -c /gpfs/hpchome/a72094/datasets/processed/GENCORD/featureCounts/merged_gene_counts.txt\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -o /gpfs/hpchome/a72094/datasets/processed/expression_matrices/