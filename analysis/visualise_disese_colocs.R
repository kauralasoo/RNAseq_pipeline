extracted_results = read.table("results/extracted_variants/results.txt", header = T, stringsAsFactors = F) %>%
  dplyr::as_tibble() %>%
  dplyr::select(qtl_group, phenotype_id, snp_id, p_nominal, beta, z_score, se) %>%
  dplyr::mutate(qtl_group = stringr::str_replace(qtl_group, "results/qtl_summary_stats/", "")) %>%
  tidyr::separate(qtl_group, into = c("study", "quant", "qtl_group"), sep = "\\/")

slfn5 = dplyr::filter(extracted_results, phenotype_id == "ENSG00000166750")

ggplot(slfn5, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()
