
#### Microarray datasets ####
#Read all metadata
fairfax_2014 = readr::read_tsv("metadata/cleaned/Fairfax_2014.tsv")
fairfax_2012 = readr::read_tsv("metadata/cleaned/Fairfax_2012.tsv")
kasela_2017 = readr::read_tsv("metadata/cleaned/Kasela_2017.tsv")
naranbhai_2015 = readr::read_tsv("metadata/cleaned/Naranbhai_2015.tsv")
CEDAR = readr::read_tsv("metadata/cleaned/CEDAR.tsv") %>%
  dplyr::mutate(batch = as.character(batch))

#Count samples in QTL groups
joint_metadata = dplyr::bind_rows(fairfax_2012, fairfax_2014, kasela_2017, naranbhai_2015, CEDAR) %>%
  dplyr::mutate(qtl_group = ifelse(qtl_group == "monocyte_CD14", "monocyte_naive", qtl_group)) %>%
  dplyr::mutate(qtl_group = ifelse(qtl_group == "neutrophil_CD16", "neutrophil", qtl_group)) %>%
  dplyr::mutate(qtl_group = ifelse(qtl_group == "neutrophil_CD15", "neutrophil", qtl_group)) %>%
  dplyr::filter(rna_qc_passed, genotype_qc_passed)

#Count samples per group
counts = table(joint_metadata$qtl_group)
counts_df = data_frame(qtl_group = names(counts), count = counts)
write.table(counts_df,"results/QC/HumanHT-12_V4.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


#### RNA-seq datasets ####
GENCORD = readr::read_tsv("metadata/cleaned/GENCORD.tsv")
BLUEPRINT = readr::read_tsv("metadata/cleaned/BLUEPRINT.tsv")
GEUVADIS = readr::read_tsv("metadata/cleaned/GEUVADIS.tsv")
TwinsUK = readr::read_tsv("metadata/cleaned/TwinsUK.tsv")
Quach_2016 = readr::read_tsv("metadata/cleaned/Quach_2016.tsv")
Nedelec_2016 = readr::read_tsv("metadata/cleaned/Nedelec_2016.tsv")

joint_metadata = dplyr::bind_rows(GENCORD, BLUEPRINT, GEUVADIS, TwinsUK, Quach_2016, Nedelec_2016) %>%
  dplyr::filter(rna_qc_passed, genotype_qc_passed) %>%
  dplyr::mutate(qtl_group = paste(cell_type, condition, sep = "_"))
table(joint_metadata$qtl_group)

