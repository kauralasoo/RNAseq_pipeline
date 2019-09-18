library("dplyr")
library("purrr")

importColocTable <- function(coloc_dir) {
  coloc_path = file.path(coloc_dir)
  full_paths = list.files(coloc_path, full.names = T)
  short_paths = list.files(coloc_path)
  
  #Extract qtl_group
  paths_df = dplyr::data_frame(file_name = full_paths, short_name = short_paths) %>%
    tidyr::separate(short_name, c("study","qtl_group", "quant", "gwas_trait", "suffix"), sep = "\\.") %>%
    dplyr::select(study, qtl_group, quant, gwas_trait, file_name) %>%
    dplyr::distinct() 

  #All files
  all_files = setNames(full_paths, full_paths)
  imported_files = purrr::map_df(all_files, ~read.table(., header = T, stringsAsFactors = F) %>% dplyr::as_tibble(), .id = "file_name") %>%
    dplyr::left_join(paths_df, by = c("file_name", "gwas_trait")) %>%
    dplyr::select(-file_name) %>%
    dplyr::mutate(PP_power = PP.H3.abf + PP.H4.abf, PP_coloc = PP.H4.abf/(PP.H3.abf + PP.H4.abf)) 
  return(imported_files)
}

#Import metadata
array_meta = read.table("~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz", header = T, stringsAsFactors = F) %>%
  dplyr::select(phenotype_id, gene_id, gene_name)

#Import Astle gwas names
astle_gwas_names = read.table("metadata/gwas_metadata/Astle_2016.csv", sep = ",", header = T, stringsAsFactors = F)[,c("Study.accession", "Trait.s.")]
colnames(astle_gwas_names) = c("trait_id", "trait_name")

#Array coloc results
folders = c("results/coloc_v2/CEDAR/", "results/coloc_v2/Fairfax_2014/", "results/coloc_v2/Fairfax_2012/", 
            "results/coloc_v2/Naranbhai_2015/", "results/coloc_v2/Kasela_2017/")
coloc_results = purrr::map_df(setNames(as.list(folders), folders), ~importColocTable(.)) %>% 
  dplyr::left_join(array_meta, by = "phenotype_id") %>%
  tidyr::separate(gwas_trait, c("pubmed_id", "trait_id", "efo_id"), "-", keep = T) %>%
  dplyr::left_join(astle_gwas_names, by = "trait_id")

hits = dplyr::filter(coloc_results, PP_power > 0.8, PP_coloc > 0.9) %>% 
  dplyr::select(gene_name, PP_power, PP_coloc, trait_name, qtl_group)
write.table(hits, "results/coloc_v2/coloc_hits.txt", sep = "\t", quote = F, row.names = F)


