
sumstats_home = "/gpfs/hpc/home/a72094/datasets/summary_stats/eQTLCatalogue/v0.1/"
matrix_home = "/gpfs/hpc/home/a72094/datasets/processed/expression_matrices/normalised_studies/"
sample_meta_home = "/gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/"
phenotype_meta_home = "/gpfs/hpc/home/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/"

quants = c("ge","exon","tx", "txrev")
quant_method = c("gene_counts", "exon_counts", "transcript_usage", "txrevise")
matrix_suffix = c("gene_counts_cqn_norm.tsv","exon_counts_cqn_norm.tsv","transcript_usage_qnorm.tsv", "txrevise_qnorm.tsv")
meta_file_name = c("gene_counts_Ensembl_96_phenotype_metadata.tsv.gz","exon_counts_Ensembl_96_phenotype_metadata.tsv.gz","transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz","txrevise_Ensembl_96_phenotype_metadata.tsv.gz")

#Study speciic
study = "Alasoo_2018"
qtl_groups = c("macrophage_naive", "macrophage_IFNg", "macrophage_Salmonella", "macrophage_IFNg+Salmonella")
vcf = "/gpfs/hpc/home/a72094/datasets/controlled_access/Alasoo_2018/genotypes/Alasoo_2018_GRCh38.filtered.vcf.gz"

qtl = dplyr::tibble(qtl_group = qtl_groups, c = 1)
quant = dplyr::tibble(quant_method = quant_method, quant = quants, matrix_suffix = matrix_suffix, meta_file_name = meta_file_name, c = 1)
susie_file = dplyr::left_join(qtl,quant, by = "c") %>%
  dplyr::select(-c) %>%
  dplyr::mutate(study = study) %>%
  dplyr::mutate(set = paste(study, quant, qtl_group, sep = "_")) %>%
  dplyr::mutate(covariates = file.path(sumstats_home, "pipeline_output", study, "PCA", set, paste0(set, ".covariates.txt"))) %>%
  dplyr::mutate(phenotype_list = file.path(sumstats_home, "summary_stats", study,paste0(set,".permuted.txt.gz"))) %>%
  dplyr::mutate(expression_matrix = file.path(matrix_home, study, paste0(study, ".", matrix_suffix))) %>%
  dplyr::mutate(sample_meta = file.path(sample_meta_home, paste0(study, ".tsv"))) %>%
  dplyr::mutate(vcf = vcf) %>%
  dplyr::mutate(phenotype_meta = file.path(phenotype_meta_home, meta_file_name)) %>%
  dplyr::select(study, qtl_group, quant_method, expression_matrix, phenotype_meta, sample_meta, vcf, phenotype_list, covariates)

susie_file = dplyr::filter(susie_file, quant_method != "exon_counts")
write.table(susie_file, "../SampleArcheology/finemapping/Alasoo_2018.tsv", sep = "\t", quote = F, row.names = F)
