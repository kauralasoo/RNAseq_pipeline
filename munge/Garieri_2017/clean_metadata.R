library("dplyr")
library("stringr")

#Import GENCORD metadata
GENCORD_meta = read.table("metadata/cleaned/GENCORD.tsv", header = TRUE, stringsAsFactors = FALSE) %>%
  dplyr::transmute(new_gt = genotype_id, genotype_id = sample_id)

meta = read.table("metadata/Garieri_2017/PRJEB21597.txt", stringsAsFactors = FALSE, header = TRUE, sep = "\t") %>%
  as_tibble() %>%
  dplyr::transmute(sample_id = run_accession, ENA_sample_title = sample_title) %>%
  tidyr::separate(ENA_sample_title, c("s1", "s3", "sequence_tag", "genotype_id"), sep = "_", remove = F) %>%
  dplyr::mutate(batch_id = paste(s1, s3, sep = "_")) %>%
  dplyr::select(sample_id, genotype_id, batch_id, sequence_tag, ENA_sample_title) %>%
  dplyr::left_join(GENCORD_meta, by = "genotype_id") %>%
  dplyr::mutate(new_gt = ifelse(is.na(new_gt), genotype_id, new_gt)) %>%
  dplyr::mutate(genotype_id = new_gt) %>%
  dplyr::select(-new_gt)

write.table(meta, "metadata/Garieri_2017/Garieri_sample_metadata.txt", sep = "\t", quote = FALSE, row.names = F)
write.table(meta$genotype_id, "metadata/Garieri_2017/Garieri_genotype_names.txt", sep = "\t", quote = FALSE, row.names = F, col.names = F)
