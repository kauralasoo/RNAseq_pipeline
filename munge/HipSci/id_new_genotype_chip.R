library("dplyr")

#Find high missingness
missingness = read.table("metadata/HipSci/HipSci_12v1-0_GRCh37_missing.imiss", header = T, stringsAsFactors = F)
missing_indiv = dplyr::filter(missingness, F_MISS > 0.3) %>%
  dplyr::rename(Sample_ID = INDV)

#Load GenomeStudio config file
config = readr::read_delim("metadata/HipSci/GenomeStudio_sample_sheet.csv", skip = 8, delim = ",")

#Extract new chip
new_chip = dplyr::semi_join(config, missing_indiv, by = "Sample_ID")
write.table(new_chip, "metadata/HipSci/GenomeStudio_HumanCoreExome12-v1-1-C.csv", sep = ",", quote = F, row.names = F)

#Extract old chip
old_chip = dplyr::anti_join(config, missing_indiv, by = "Sample_ID")
write.table(old_chip, "metadata/HipSci/GenomeStudio_HumanCoreExome12-v1-0-D.csv", sep = ",", quote = F, row.names = F)
