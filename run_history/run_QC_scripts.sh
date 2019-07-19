#GENCORD
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpchome/a72094/datasets/processed/GENCORD/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir ~/datasets/processed/GENCORD/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/GENCORD/

#Schmiedel_2018
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpchome/a72094/projects/Schmiedel_2018/rnaseq/results/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Schmiedel_2018.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpchome/a72094/projects/Schmiedel_2018/rnaseq/results/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Schmiedel_2018/

#TwinsUK
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpchome/a72094/projects/rnaseq/results/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/TwinsUK.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpchome/a72094/projects/rnaseq/results/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/TwinsUK/\
  --build_html TRUE

#Lepik_2017
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpchome/a72094/datasets/processed/Lepik_2017/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Lepik_2017.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpchome/a72094/datasets/processed/Lepik_2017/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Lepik_2017/\
  --build_html TRUE
