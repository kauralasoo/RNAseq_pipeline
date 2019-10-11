library("dplyr")
library("purrr")

coloc_res = read.table("results/full_coloc_res.tsv", sep = "\t", header = TRUE, stringsAsFactors = F) %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(PP_power = PP.H3.abf + PP.H4.abf, PP_coloc = PP.H4.abf/(PP.H3.abf + PP.H4.abf)) %>%
  dplyr::filter(PP.H4.abf > 0.5) %>%
  #dplyr::filter(PP_power > 0.8, PP_coloc > 0.9) %>%
  dplyr::group_by(snp_id, gwas_trait_id) %>%
  dplyr::arrange(snp_id,gwas_trait_id, -PP.H4.abf) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()
View(coloc_res)

write.table(coloc_res, "results/coloc_results.txt", sep = "\t", quote = F, row.names = F)

a = read.table("~/Downloads/full_cs_matchdata.tsv", header = T, sep = "\t", stringsAsFactors = F)

hits = dplyr::select(a, meta_id, cell_type, study, cis_phenotype, gene_name) %>% dplyr::distinct()
