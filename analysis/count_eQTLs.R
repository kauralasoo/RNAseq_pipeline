library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("ggplot2")
load_all("../eQTLUtils/")
load_all("../seqUtils/")

#Import all lead variants
lead_vars = idVectorToList(list.files("results/summary_stats/", recursive = T, full.names = T))
lead_var_list = purrr::map(lead_vars, ~eQTLUtils::importQTLtoolsTable(.))

#Count QTLs
qtl_count = purrr::map_df(lead_var_list, ~dplyr::filter(., p_fdr < 0.05) %>% dplyr::summarise(eQTL_count = length(phenotype_id)), .id = "path")
feature_count = purrr::map_df(lead_var_list, ~dplyr::summarise(.,feature_count = length(phenotype_id)), .id = "path")
snp_count = purrr::map_df(lead_var_list, ~dplyr::summarise(.,mean_snp_count = mean(n_cis_snps), sd_snp_count = sd(n_cis_snps)), .id = "path")

qtl_counts = dplyr::left_join(qtl_count, feature_count, by = "path") %>% 
  dplyr::left_join(snp_count, by = "path") %>%
  dplyr::mutate(path = stringr::str_replace(path,"results/summary_stats//", "")) %>%
  dplyr::mutate(path = stringr::str_replace(path, ".permuted.txt.gz", "")) %>%
  tidyr::separate(path, c("study", "quant_method", "qtl_group"), "\\/")

#Import study metadata
meta_list = list(BLUEPRINT = "metadata/cleaned/BLUEPRINT.tsv",
                 CEDAR = "metadata/cleaned/CEDAR.tsv",
                 Fairfax_2014 = "metadata/cleaned/Fairfax_2014.tsv",
                 GENCORD = "metadata/cleaned/GENCORD.tsv",
                 GEUVADIS = "metadata/cleaned/GEUVADIS.tsv")
meta_imported = purrr::map(meta_list, ~read.table(., sep = "\t", stringsAsFactors = F, header =T) %>% dplyr::as_tibble())
samples = purrr::map_df(meta_imported, ~dplyr::select(.,cell_type, condition, qtl_group, rna_qc_passed, genotype_qc_passed, study))
cell_types = dplyr::select(samples, cell_type, condition, qtl_group, study) %>% dplyr::distinct()
sample_sizes = dplyr::filter(samples, rna_qc_passed, genotype_qc_passed) %>% 
  dplyr::group_by(study, qtl_group) %>% 
  dplyr::summarise(sample_size = length(qtl_group)) %>% 
  dplyr::ungroup()

#Merge counts
merged_counts = dplyr::left_join(qtl_counts, sample_sizes, by = c("study", "qtl_group")) %>% 
  dplyr::left_join(cell_types, by = c("qtl_group", "study"))
ggplot(merged_counts, aes(x = sample_size, y = eQTL_count, color = study, label = cell_type)) + 
  geom_point() +
  scale_y_continuous(limits = c(0, 7500)) +
  scale_x_continuous(limits = c(0, 460)) +
  ylab("Number of eGenes") +
  xlab("sample size")

#Visualise SNP counts per study
snp_counts_df = dplyr::select(qtl_counts, study, mean_snp_count, sd_snp_count) %>% dplyr::distinct()
ggplot(snp_counts_df, aes(x = study, y = mean_snp_count, ymax = mean_snp_count+sd_snp_count, ymin = mean_snp_count-sd_snp_count, fill = study)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(width = 0.2)

