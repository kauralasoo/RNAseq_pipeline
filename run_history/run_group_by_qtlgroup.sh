#GENCORD featureCounts
singularity exec work/singularity/kerimoff-qtlmap-latest.img Rscript bin/group_by_qtlgroup.R\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -e /gpfs/hpchome/a72094/datasets/processed/expression_matrices/GENCORD_gene_counts_cqn_norm.tsv\
 -v /gpfs/hpchome/a72094/datasets/controlled_access/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.variant_information.txt.gz\
 -t /gpfs/hpchome/a72094/datasets/processed/QC_reports/GENCORD/median_tpm/GENCORD_95quantile_tpm.tsv.gz\
 -o test/

#GENCORD leafcutter
 singularity exec work/singularity/kerimoff-qtlmap-latest.img Rscript bin/group_by_qtlgroup.R\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p qtl_data/UCT_phenotype_metadata.tsv\
 -e qtl_data/UCT_qnorm.tsv\
 -v /gpfs/hpchome/a72094/datasets/controlled_access/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.variant_information.txt.gz\
 -o test/\
 -t /gpfs/hpchome/a72094/datasets/processed/QC_reports/GENCORD/median_tpm/GENCORD_95quantile_tpm.tsv.gz