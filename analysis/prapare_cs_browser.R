library("dplyr")

#Extract effect size from sumstats
tabix_table = readr::read_tsv("results/credible_set_browser/tabix_ftp_paths.tsv") %>%
  dplyr::filter(quant_method %in% c("ge","microarray")) %>%
  dplyr::mutate(qtl_subset = paste0(study, ".", qtl_group, "_", quant_method)) %>%
  dplyr::mutate(credible_sets = paste0("/gpfs/space/home/a72094/eQTLCatalogue/susie-finemapping/eQTL_Catalogue_r3/susie/",qtl_subset, ".purity_filtered.txt.gz")) %>%
  dplyr::mutate(qtl_ss = paste0("/gpfs/space/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/",qtl_subset, ".nominal.sorted.tsv.gz")) %>%
  dplyr::mutate(qtl_ss_index = paste0("/gpfs/space/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out/sumstats/",qtl_subset, ".nominal.sorted.tsv.gz.tbi"))

extract_table = dplyr::select(tabix_table, qtl_subset, credible_sets, qtl_ss, qtl_ss_index)
write.table(extract_table, "results/credible_set_browser/effect_extract.tsv", sep = "\t", row.names = F, quote = F)

#Make the final credible set table
dataset_table = dplyr::mutate(tabix_table, dataset = paste(study, tissue_label, condition_label, sep = ":")) %>% 
  dplyr::select(qtl_subset, dataset) %>%
  dplyr::mutate(dataset_index = paste0("D", c(1:length(dataset)))) %>%
  dplyr::filter(qtl_subset != "GTEx.kidney_cortex_ge")

importCsData <- function(qtl_subset_id){
  cs_path = paste0("results/credible_set_browser/susie/", qtl_subset_id, ".purity_filtered.txt.gz")
  cs_df = readr::read_tsv(cs_path) %>% 
    dplyr::transmute(molecular_trait_id = phenotype_id, variant = variant_id, credible_set = cs_id, pip, credible_set_size = cs_size) %>%
    dplyr::mutate(qtl_subset = qtl_subset_id)
  
  sumstats_path = paste0("results/credible_set_browser/extracted_sumstats/", qtl_subset_id, ".extracted_sumstats.tsv.gz")
  sumstats_df = readr::read_tsv(sumstats_path)
  overlapping_cs = dplyr::left_join(cs_df, sumstats_df, by = c("molecular_trait_id", "variant"))
  
  return(overlapping_cs)
}

dataset_list = setNames(as.list(dataset_table$qtl_subset), dataset_table$qtl_subset)
#sumstats_table = purrr::map_df(dataset_list, importCsData)
#saveRDS(sumstats_table, "results/credible_set_browser/sumstats_table.rds")
sumstats_table = readRDS("results/credible_set_browser/sumstats_table.rds")

#Import gene and probe metadata
#Import gene metadata
gene_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::select(phenotype_id, gene_name)
probe_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz") %>%
  dplyr::select(phenotype_id, gene_name)
gene_map = dplyr::bind_rows(gene_meta, probe_meta) %>%
  dplyr::rename(molecular_trait_id = phenotype_id)

#Make the final table
final_table = dplyr::left_join(sumstats_table, dataset_table, by = "qtl_subset") %>%
  dplyr::left_join(gene_map, by = "molecular_trait_id")

save_table = dplyr::mutate(final_table, credible_set = paste0(credible_set, "_", dataset_index)) %>%
  dplyr::select(molecular_trait_id, gene_name, credible_set, variant, rsid, credible_set_size, pip, pvalue, beta, se, dataset)
file_handle = gzfile("results/credible_set_browser/credible_set_table.tsv.gz","w")
write.table(save_table, file_handle, sep = "\t", row.names = F, col.names = T, quote = FALSE)
close(file_handle)


