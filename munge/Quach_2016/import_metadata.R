sample_meta = read.table("metadata/Quach_2016/file_names.txt", stringsAsFactors = FALSE) %>%
  dplyr::tbl_df() %>%
  dplyr::rename(fq1 = V1) %>%
  tidyr::separate(fq1, c("prefix","ega_id", "donor_id", "suffix"), sep = "_", remove = FALSE) %>%
  dplyr::select(-prefix, -suffix) %>%
  tidyr::separate(donor_id, c("pre", "sample"), sep = "EIP@") %>%
  tidyr::separate(sample, c("donor", "condition_number"), sep = "-") %>%
  dplyr::mutate(condition_name = case_when(
    condition_number == 1 ~ "naive",
    condition_number == 2 ~ "LPS",
    condition_number == 3 ~ "Pam3CSK4",
    condition_number == 4 ~ "R848",
    condition_number == 5 ~ "IAV"
   )) %>%
  dplyr::mutate(sample_id = paste(donor, condition_number, sep = "_")) %>%
  dplyr::select(sample_id, donor, condition_number, condition_name, ega_id, fq1)

write.table(sample_meta, "metadata/Quach_2016/Quach_2016_compiled_metadata.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Make fastq paths
snakemake = dplyr::transmute(sample_meta, sample_id, fq1) %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/Monocytes_Quach_2016/fastq/", fq1,"]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakemake, "metadata/Quach_2016/Quach_2016_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
