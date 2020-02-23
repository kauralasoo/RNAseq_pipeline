library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
load_all("../eQTLUtils/")

#Share with Melanie
meta_list = list(Schwartzentruber_2018 = "../SampleArcheology/studies/cleaned/Schwartzentruber_2018.tsv",
                 van_de_Bunt_2015 = "../SampleArcheology/studies/cleaned/van_de_Bunt_2015.tsv",
                 HipSci = "../SampleArcheology/studies/cleaned/HipSci.tsv")
meta_table = purrr::map_df(meta_list, ~read.table(., sep = "\t", stringsAsFactors = F, header =T) %>% 
                                dplyr::as_tibble() %>%
                                dplyr::transmute(sample_id, genotype_id, cell_type, study, protocol, paired, stranded, read_length))
write.table(meta_table, "results/hipsci_bunt_schwartzentruber_sample_ids.txt", sep = "\t", quote = F, row.names = F)

#Share with Taru
#Import metadata files
meta_list = list(Schwartzentruber_2018 = "../SampleArcheology/studies/cleaned/Schwartzentruber_2018.tsv",
                 van_de_Bunt_2015 = "../SampleArcheology/studies/cleaned/van_de_Bunt_2015.tsv",
                 HipSci = "../SampleArcheology/studies/cleaned/HipSci.tsv",
                 BLUEPRINT_SE = "../SampleArcheology/studies/cleaned/BLUEPRINT_SE.tsv",
                 BLUEPRINT_PE = "../SampleArcheology/studies/cleaned/BLUEPRINT_PE.tsv",
                 GEUVADIS = "../SampleArcheology/studies/cleaned/GEUVADIS.tsv",
                 Alasoo_2018 = "../SampleArcheology/studies/cleaned/Alasoo_2018.tsv",
                 Schmiedel_2018 = "../SampleArcheology/studies/cleaned/Schmiedel_2018.tsv",
                 ROSMAP = "../SampleArcheology/studies/cleaned/ROSMAP.tsv",
                 BrainSeq = "../SampleArcheology/studies/cleaned/BrainSeq.tsv",
                 Fairfax_2014 = "../SampleArcheology/studies/cleaned/Fairfax_2014.tsv",
                 Kasela_2017 = "../SampleArcheology/studies/cleaned/Kasela_2017.tsv",
                 Fairfax_2012 = "../SampleArcheology/studies/cleaned/Fairfax_2012.tsv",
                 Naranbhai_2015 = "../SampleArcheology/studies/cleaned/Naranbhai_2015.tsv",
                 CEDAR = "../SampleArcheology/studies/cleaned/CEDAR.tsv"
)
meta_table = purrr::map_df(meta_list, ~read.table(., sep = "\t", stringsAsFactors = F, header =T) %>% 
                             dplyr::as_tibble() %>%
                             dplyr::transmute(sample_id, sex, cell_type, condition, qtl_group, study, rna_qc_passed, genotype_qc_passed, protocol, paired, stranded, read_length))
write.table(meta_table, "results/sex_differences_metadata.txt", sep = "\t", quote = F, row.names = F)
