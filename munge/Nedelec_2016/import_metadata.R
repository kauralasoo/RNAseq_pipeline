library("dplyr")

metadata = readr::read_csv("metadata/Nedelec_2017/Nedelec_2017_metadata.csv")
sample_names = readr::read_csv("metadata/Nedelec_2017/Nedelec_2017_sample_names.csv") %>%
  dplyr::rename(geo_sample_id = Sample_geo_accession, sample_id = Sample_title)

sample_meta = dplyr::transmute(metadata, run_id = Run, geo_sample_id = SampleName) %>% 
  dplyr::left_join(sample_names) %>%
  dplyr::select(sample_id, everything())

#Import metadata from the article
article_samples = readr::read_csv("metadata/Nedelec_2017/article_sample_metadata.csv") %>%
  tidyr::separate(ID, c("sample_id", "NN"), sep = "_mRNA") %>%
  dplyr::select(-NN) %>%
  dplyr::rename(donor = Label, condition = Condition, sequencing_center = Sequencing_center, 
                flowcell = Flowcell, lane = Lane) %>%
  dplyr::mutate(condition_name = case_when(
    condition == "Listeria" ~ "L",
    condition == "Non-infected" ~ "NI",
    condition == "Salmonella" ~ "S"
  ))

article_donors = readr::read_csv("metadata/Nedelec_2017/article_donor_metadata.csv") %>%
  dplyr::transmute(donor = Label, self_reported_ethnic = Self_reported_ethnic, AF_admixture = African_admixture,
                   EU_admixture = European_admixture, ethnic_bin = Ethnic_bin)

joint_meta = dplyr::left_join(sample_meta, article_samples, by = "sample_id") %>% 
  dplyr::left_join(article_donors, by = "donor")
write.table(joint_meta, "metadata/Nedelec_2017/Nedelec_2016_compiled_metadata.txt", 
            sep = "\t", quote = FALSE, row.names = FALSE)

#Make fastq paths
snakemake = dplyr::transmute(joint_meta, sample_id, fq1 = run_id) %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/Macrophages_Nedelec_2016/fastq/", fq1,"_1.fastq.gz]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakemake, "metadata/Nedelec_2017/Nedelec_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

