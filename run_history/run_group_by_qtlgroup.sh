#GENCORD featureCounts
singularity exec work/singularity/kerimoff-qtlmap-latest.img Rscript bin/group_by_qtlgroup.R\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -e /gpfs/hpchome/a72094/datasets/processed/expression_matrices/GENCORD_gene_counts_cqn_norm.tsv\
 -v /gpfs/hpchome/a72094/datasets/controlled_access/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.variant_information.txt.gz\
 -t /gpfs/hpchome/a72094/datasets/processed/QC_reports/GENCORD/median_tpm/GENCORD_95quantile_tpm.tsv.gz\
 -o test/

 #GENCORD transcript_usage
singularity exec work/singularity/kerimoff-qtlmap-latest.img Rscript bin/group_by_qtlgroup.R\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/GENCORD.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz\
 -e /gpfs/hpchome/a72094/datasets/processed/transcript_usage/GENCORD.transcript_usage_qnorm.tsv\
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


 #Lepik_2017 exon_counts
singularity exec work/singularity/kerimoff-qtlmap-latest.img Rscript bin/group_by_qtlgroup.R\
 -s /gpfs/hpchome/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/Lepik_2017.tsv\
 -p /gpfs/hpchome/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/exon_counts_Ensembl_96_phenotype_metadata.tsv.gz\
 -e /gpfs/hpchome/a72094/datasets/processed/expression_matrices/exon_counts/Lepik_2017_exon_counts_cqn_norm.tsv\
 -v /gpfs/hpchome/a72094/datasets/controlled_access/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.variant_information.txt.gz\
 -t /gpfs/hpchome/a72094/datasets/processed/QC_reports/Lepik_2017/median_tpm/Lepik_2017_95quantile_tpm.tsv.gz\
 -o test/
 