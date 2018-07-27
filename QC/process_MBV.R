library("stringr")
library("ggplot2")

#Get file names
mbv_files = list.files("processed/GEUVADIS/mbv/", full.names = T)

#Make sample names
sample_names = stringr::str_replace_all(basename(mbv_files), ".mbv_output.txt", "")
sample_list = setNames(mbv_files, sample_names)

#Import mbv files
mbv_results = purrr::map(sample_list, ~readr::read_delim(., delim = " ", col_types = "ciiiiiiiiii"))

findBestMatch <- function(mbv_df){
  res = dplyr::transmute(mbv_df, mbv_genotype_id = SampleID, 
                         het_consistent_frac = n_het_consistent/n_het_covered, 
                         hom_consistent_frac = n_hom_consistent/n_hom_covered)
  
  #Identify best het
  best_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  other_het = dplyr::arrange(res, -het_consistent_frac) %>% dplyr::filter(dplyr::row_number() > 1)
  best_row = dplyr::mutate(best_het, het_min_dist = min(best_het$het_consistent_frac - other_het$het_consistent_frac),
                           hom_min_dist = min(best_het$hom_consistent_frac - other_het$hom_consistent_frac))
  
  #Compare against best hom
  best_hom = dplyr::arrange(res, -hom_consistent_frac) %>% dplyr::filter(dplyr::row_number() == 1)
  if(best_row$mbv_genotype_id != best_hom$mbv_genotype_id){
    best_row = dplyr::mutate(best_row, het_consistent_frac = as.numeric(NA), hom_consistent_frac = as.numeric(NA),
                  het_min_dist = as.numeric(NA), hom_min_dist = as.numeric(NA))
  }
  return(best_row)
}

best_matches = purrr::map_df(mbv_results, findBestMatch, .id = "sample_id") %>%
  dplyr::filter(!is.na(het_consistent_frac)) %>%
  dplyr::filter(het_consistent_frac > 0.9)
write.table(best_matches, "metadata/GEUVADIS/GEUVADIS_mbv_best_match.txt", sep = "\t", quote = F, row.names = F)


