library("dplyr")

#Import raw files
sra_run_table = read.table("metadata/Ye_2018/SraRunTable.txt", sep = "\t", header = T, stringsAsFactors = F) %>%
  dplyr::as_tibble()
srr_samples = read.table("metadata/Ye_2018/ImmVar_SRA_samples.txt", skip = 10, header = T, sep = "\t")

#Extract dendritic cells
DC_samples = dplyr::filter(sra_run_table, Sample_Name %in% srr_samples$SAMPLE_ID) %>%
  dplyr::transmute(biosample_id = BioSample, sample_id = Run, sample_name = Sample_Name, sex, donor_id = submitted_subject_id)
write.table(DC_samples, "metadata/Ye_2018/Ye_2018_dencritic_cells.txt", sep = "\t", quote = F, row.names = F)

