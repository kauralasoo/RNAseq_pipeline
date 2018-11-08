library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cqn")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
library("devtools")
load_all("../eQTLUtils/")

#Specify mandatory metadata columns
mandatory_cols = c("sample_id", "genotype_id", "sex", "cell_type", "condition", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")

#Fairfax_2014
#Import SE
fairfax_2014 = readRDS("results/SummarizedExperiments/microarray/Fairfax_2014.rds")

#Import Plink sex estimates
set1 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set1.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
set2 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set2.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
plink_sex = dplyr::bind_rows(set1, set2) %>% 
  dplyr::select(genotype_id, plink_sex) %>% 
  dplyr::mutate(sex = "NA") %>% 
  dplyr::mutate(sex = ifelse(plink_sex == 1, "male", sex)) %>%
  dplyr::mutate(sex = ifelse(plink_sex == 2, "female", sex)) %>%
  dplyr::select(genotype_id, sex)

#Reformat 
sample_metadata = colData(fairfax_2014) %>% as.data.frame() %>% as_tibble() %>%
  dplyr::select(-sex) %>%
  dplyr::left_join(plink_sex, by = "genotype_id") %>%
  dplyr::rename(genotype_qc_passed = genotype_QC_passed) %>%
  dplyr::rename(rna_qc_passed = RNA_QC_passed) %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(sex), FALSE, TRUE)) %>%
  dplyr::mutate(genotype_id = paste0("FF14_", genotype_id)) %>%
  dplyr::mutate(condition = ifelse(condition == "naive", condition, paste0(condition, timepoint))) %>%
  dplyr::mutate(cell_type = "monocyte") %>%
  dplyr::mutate(protocol = "HumanHT-12_V4") %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/Fairfax_2014.tsv", sep = "\t", quote = FALSE, row.names = FALSE)


