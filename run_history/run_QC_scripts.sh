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
  -c /gpfs/hpc/home/a72094/datasets/processed/Schmiedel_2018/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Schmiedel_2018.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/Schmiedel_2018//MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Schmiedel_2018/\
  --build_html TRUE

#TwinsUK
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/TwinsUK/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/TwinsUK.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/TwinsUK/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/TwinsUK/\
  --build_html TRUE

#Nedelec_2016
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/Nedelec_2016/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Nedelec_2016.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/Nedelec_2016/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Nedelec_2016/\
  --build_html TRUE

#Quach_2016
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/Quach_2016/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Quach_2016.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/Quach_2016/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Quach_2016/\
  --build_html TRUE

#GEUVADIS
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/GEUVADIS/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GEUVADIS.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/GEUVADIS/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/GEUVADIS/\
  --build_html TRUE

#HipSci
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/HipSci/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/HipSci.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/HipSci/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/HipSci/\
  --build_html TRUE

#Schwartzentruber_2018
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/Schwartzentruber_2018/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Schwartzentruber_2018.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/Schwartzentruber_2018/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Schwartzentruber_2018/\
  --build_html TRUE

#van_de_Bunt_2015
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/van_de_Bunt_2015/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/van_de_Bunt_2015.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/van_de_Bunt_2015/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/van_de_Bunt_2015/\
  --build_html TRUE

#Lepik_2017
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/Lepik_2017/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Lepik_2017.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/Lepik_2017/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Lepik_2017/\
  --build_html TRUE

#BLUEPRINT_SE
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/BLUEPRINT_SE/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/BLUEPRINT.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/BLUEPRINT_SE/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/BLUEPRINT_SE/\
  --build_html TRUE

#BLUEPRINT_PE
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/BLUEPRINT_PE/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/BLUEPRINT.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/BLUEPRINT_PE/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/BLUEPRINT_PE/\
  --build_html TRUE


#Alasoo_2018
singularity exec qtl_norm_qc.img Rscript feature_counts_qc.R\
  -c /gpfs/hpc/home/a72094/datasets/processed/Alasoo_2018/featureCounts/merged_gene_counts.txt\
  -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
  -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
  -q ge\
  --mbvdir /gpfs/hpc/home/a72094/datasets/processed/Alasoo_2018/MBV/\
  --filter_qc TRUE\
  -o /gpfs/hpchome/a72094/datasets/processed/QC_reports/Alasoo_2018/\
  --build_html TRUE