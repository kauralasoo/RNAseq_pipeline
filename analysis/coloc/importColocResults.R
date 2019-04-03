library("dplyr")

importColocTable <- function(coloc_dir, study_id, quant_id) {
  coloc_path = file.path(coloc_dir, study_id, quant_id)
  full_paths = list.files(coloc_path, full.names = T)
  short_paths = list.files(coloc_path)
  
  #Extract qtl_group
  paths_df = dplyr::data_frame(file_name = full_paths, short_name = short_paths) %>%
    tidyr::separate(short_name, c("qtl_group", "gwas_trait", "suffix"), sep = "\\.") %>%
    dplyr::select(file_name, qtl_group) %>%
    dplyr::distinct() %>%
    dplyr::mutate(study = study_id, quant = quant_id)
  
  #All files
  all_files = setNames(full_paths, full_paths)
  imported_files = purrr::map_df(all_files, ~read.table(., header = T, stringsAsFactors = F) %>% dplyr::as_tibble(), .id = "file_name") %>%
    dplyr::left_join(paths_df, by = "file_name") %>%
    dplyr::select(-file_name, -.row) %>%
    dplyr::mutate(PP_power = PP.H3.abf + PP.H4.abf, PP_coloc = PP.H4.abf/(PP.H3.abf + PP.H4.abf)) 
  return(imported_files)
}

#Import GWAS traits
gwas_stats_labeled = readr::read_tsv("metadata/gwas_metadata/GWAS_summary_stat_list.labeled.txt",
                                     col_names = c("trait","file_name", "type")) %>%
  dplyr::filter(!(trait %in% c("UC_2014","UC_2012", "CEL_2010","PS", "CD_2012", "RA_2012", "T2D_1", "MS", "T1D", "T1D_2", "PBC")))

#Import gene metadata
array_gene_meta = read.table("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz", header = TRUE, stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
probe_gene_map = dplyr::select(array_gene_meta, phenotype_id, gene_name)
featureCounts_gene_meta = read.table("metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz", header = TRUE, stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
gene_name_map = dplyr::select(featureCounts_gene_meta, phenotype_id, gene_name)
leafcutter_gene_meta = read.table("metadata/gene_metadata/GEUVADIS_leafcutter_cluster_metadata.txt.gz", header = TRUE, stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
leafcutter_name_map = dplyr::select(leafcutter_gene_meta, phenotype_id, gene_name)


#Import array QTL results
array_studies = c("Fairfax_2012", "Fairfax_2014", "Kasela_2017","CEDAR","Naranbhai_2015")
array_list = setNames(array_studies, array_studies)

array_colocs = purrr::map_df(array_list, ~importColocTable("results/coloc/", ., "array")) %>%
  dplyr::left_join(probe_gene_map, by = "phenotype_id")

#Find hits
array_coloc_hits = dplyr::filter(array_colocs, PP_power > 0.8, PP_coloc > 0.9)


#Import RNA-seq studies
rnaseq_studies = c("Alasoo_2018", "Schwartzentruber_2018","BLUEPRINT","GENCORD","TwinsUK","GEUVADIS","van_de_Bunt_2015","Nedelec_2016", "HipSci", "Quach_2016")
rnaseq_list = setNames(rnaseq_studies, rnaseq_studies)

rnaseq_colocs = purrr::map_df(rnaseq_list, ~importColocTable("results/coloc/", ., "featureCounts")) %>%
  dplyr::left_join(gene_name_map, by = "phenotype_id")

#Find hits
rnaseq_coloc_hits = dplyr::filter(rnaseq_colocs, PP_power > 0.8, PP_coloc > 0.9)


#Leafcutter coloc hits
leafcutter_studies = c("GEUVADIS_EUR_100kb_coloc")
lc_list = setNames(leafcutter_studies, leafcutter_studies)

leafcutter_colocs = purrr::map_df(lc_list, ~importColocTable("results/coloc/", ., "leafcutter")) %>%
  dplyr::left_join(leafcutter_name_map, by = "phenotype_id")

#Find hits
leafcutter_coloc_hits = dplyr::filter(leafcutter_colocs, PP_power > 0.8, PP_coloc > 0.9)

#Merge all coloc together
joint_colocs = dplyr::bind_rows(array_colocs, rnaseq_colocs, leafcutter_colocs)
write.table(joint_colocs, "results/coloc/eQTLCatalogue_coloc_results.txt", sep = "\t", quote = F, row.names = F)

