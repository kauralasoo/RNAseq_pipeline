library("dplyr")
library("data.table")

#Import GENCORD metadata
ega_data = read.table("metadata/GENCORD/Sample_File.map") %>%
  dplyr::as_tibble() %>%
  dplyr::select(V1, V3, V4) %>%
  dplyr::rename(sample_id = V1, file_name = V3, ega_file_name = V4) %>%
  dplyr::mutate(cell_type = case_when(
    sample_id %like% "UCB" ~ "LCL",
    sample_id %like% "UCF" ~ "fibroblast",
    sample_id %like% "UCT" ~ "T-cell"
  )) %>%
  dplyr::mutate(condition = "naive", read_length = "50bp", 
                paired = TRUE, stranded = FALSE, 
                study = "GENCORD", protocol = "poly(A)") %>%
  dplyr::mutate(genotype_id = stringr::str_replace_all(sample_id, "UCB", "UC")) %>%
  dplyr::mutate(genotype_id = stringr::str_replace_all(genotype_id, "UCT", "UC")) %>%
  dplyr::mutate(genotype_id = stringr::str_replace_all(genotype_id, "UCF", "UC")) %>%
  dplyr::mutate(file_name = stringr::str_remove(file_name, "\\.cip"))





#### Identify missing files and re-download them ####
#Identify missing files
dled = read.table("metadata/GENCORD/GENCORD_files.txt", stringsAsFactors = FALSE) %>%
  dplyr::rename(bam_name = V1) %>%
  dplyr::mutate(file_name = stringr::str_replace(bam_name, "_EGAR\\d+_", "")) %>%
  as_tibble()

merged_data = dplyr::left_join(ega_data, dled, by = "file_name")
write.table(merged_data, "metadata/GENCORD/GENCORD_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

  
missing = dplyr::anti_join(ega_data, dled, by = "file_name") %>% 
  dplyr::select(ega_file_name)
write.table(missing, "metadata/GENCORD/missing_files.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


