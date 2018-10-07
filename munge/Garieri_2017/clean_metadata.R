library("dplyr")

meta = read.table("metadata/Garieri_2017/PRJEB21597.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::transmute(sample_id = run_accession, ENA_sample_title = sample_title) %>%
  tidyr::separate(ENA_sample_title, c("s1", "s3", "sequence_tag", "genotype_id"), sep = "_", remove = F) %>%
  dplyr::mutate(batch_id = paste(s1, s3, sep = "_")) %>%
  dplyr::select(sample_id, genotype_id, batch_id, sequence_tag, ENA_sample_title)
write.table(meta, "metadata/Garieri_2017/Garieri_sample_metadata.txt", sep = "\t", quote = FALSE, row.names = F)
