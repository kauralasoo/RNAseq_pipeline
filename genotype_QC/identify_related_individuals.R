
#### CEDAR ####
#Identify related individuals
rel_mat = read.table("~/Downloads/relatedness_matrix.tsv", header = T, stringsAsFactors = F)
rel_df = tidyr::pivot_longer(rel_mat, -iid, names_to = "iid2", values_to = "relatedness")

related_samples = dplyr::filter(rel_df, relatedness > 0.025, iid != iid2) %>% 
  dplyr::arrange(relatedness)

cedar_meta = read.table("~/projects/SampleArcheology/studies/cleaned/CEDAR.tsv", header = TRUE, stringsAsFactors = F, sep = "\t") %>% dplyr::as_tibble()

cedar_related = dplyr::filter(cedar_meta, genotype_id %in% related_samples$iid) %>%
  dplyr::arrange(genotype_qc_passed)
View(cedar_related)

blacklist = c("IPC368","IPC421","IPC045","IPC174","IPC162","IPC276","IPC103","IPC112","IPC150")
cedar_unrelated = dplyr::mutate(cedar_meta, genotype_qc_passed = ifelse(genotype_id %in% blacklist, FALSE, genotype_qc_passed))
write.table(cedar_unrelated, "~/projects/SampleArcheology/studies/cleaned/CEDAR.tsv", sep = "\t", quote = FALSE, row.names = F)



#### Kasela_2017 ####
#Identify related individuals
rel_mat = read.table("~/Downloads/Kasela_QC/relatedness_matrix.tsv", header = T, stringsAsFactors = F)
rel_df = tidyr::pivot_longer(rel_mat, -iid, names_to = "iid2", values_to = "relatedness")

related_samples = dplyr::filter(rel_df, relatedness > 0.1, iid != iid2) %>% 
  dplyr::arrange(relatedness)

cedar_meta = read.table("~/projects/SampleArcheology/studies/cleaned/Kasela_2017.tsv", header = TRUE, stringsAsFactors = F, sep = "\t") %>% dplyr::as_tibble()

cedar_related = dplyr::filter(cedar_meta, genotype_id %in% related_samples$iid) %>%
  dplyr::arrange(genotype_qc_passed)
View(cedar_related)

blacklist = c("V18228","V23612")
cedar_unrelated = dplyr::mutate(cedar_meta, genotype_qc_passed = ifelse(genotype_id %in% blacklist, FALSE, genotype_qc_passed))
write.table(cedar_unrelated, "~/projects/SampleArcheology/studies/cleaned/Kasela_2017.tsv", sep = "\t", quote = FALSE, row.names = F)


#### Fairfax_2014 ####
#Identify related individuals
rel_mat = read.table("~/Downloads/Fairfax_QC/relatedness_matrix.tsv", header = T, stringsAsFactors = F, check.names = F)
rel_df = tidyr::pivot_longer(rel_mat, -iid, names_to = "iid2", values_to = "relatedness")

related_samples = dplyr::filter(rel_df, relatedness > 0.1, iid != iid2) %>% 
  dplyr::arrange(relatedness)

cedar_meta = read.table("~/projects/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv", header = TRUE, stringsAsFactors = F, sep = "\t") %>% dplyr::as_tibble()

cedar_related = dplyr::filter(cedar_meta, genotype_id %in% related_samples$iid) %>%
  dplyr::arrange(genotype_qc_passed)
View(cedar_related)

blacklist = c("FF14_16")
cedar_unrelated = dplyr::mutate(cedar_meta, genotype_qc_passed = ifelse(genotype_id %in% blacklist, FALSE, genotype_qc_passed))
write.table(cedar_unrelated, "~/projects/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv", sep = "\t", quote = FALSE, row.names = F)
