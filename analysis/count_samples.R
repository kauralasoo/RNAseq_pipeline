library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
load_all("../eQTLUtils/")

#Import metadata files
meta_list = list(BLUEPRINT = "../SampleArcheology/studies/cleaned/BLUEPRINT.tsv",
                 CEDAR = "../SampleArcheology/studies/cleaned/CEDAR.tsv",
                 Fairfax_2014 = "../SampleArcheology/studies/cleaned/Fairfax_2014.tsv",
                 GENCORD = "../SampleArcheology/studies/cleaned/GENCORD.tsv",
                 GEUVADIS = "../SampleArcheology/studies/cleaned/GEUVADIS.tsv",
                 TwinsUK = "../SampleArcheology/studies/cleaned/TwinsUK.tsv",
                 Alasoo_2018 = "../SampleArcheology/studies/cleaned/Alasoo_2018.tsv",
                 Nedelec_2016 = "../SampleArcheology/studies/cleaned/Nedelec_2016.tsv",
                 Quach_2016 = "../SampleArcheology/studies/cleaned/Quach_2016.tsv",
                 Fairfax_2012 = "../SampleArcheology/studies/cleaned/Fairfax_2012.tsv",
                 Schwartzentruber_2018 = "../SampleArcheology/studies/cleaned/Schwartzentruber_2018.tsv",
                 van_de_Bunt_2015 = "../SampleArcheology/studies/cleaned/van_de_Bunt_2015.tsv",
                 HipSci = "../SampleArcheology/studies/cleaned/HipSci.tsv",
                 Naranbhai_2015 = "../SampleArcheology/studies/cleaned/Naranbhai_2015.tsv",
                 Kasela_2017 = "../SampleArcheology/studies/cleaned/Kasela_2017.tsv",
                 Ye_2018 = "../SampleArcheology/studies/cleaned/Ye_2018.tsv",
                 Raj_2014 = "../SampleArcheology/studies/cleaned/Raj_2014.tsv",
                 Schmiedel_2018 = "../SampleArcheology/studies/cleaned/Schmiedel_2018.tsv",
                 ROSMAP = "../SampleArcheology/studies/cleaned/ROSMAP.tsv",
                 BrainSeq = "../SampleArcheology/studies/cleaned/BrainSeq.tsv")
meta_imported = purrr::map_df(meta_list, ~read.table(., sep = "\t", stringsAsFactors = F, header =T) %>% 
                             dplyr::as_tibble() %>%
                             dplyr::select(sample_id, genotype_id, cell_type, condition, qtl_group, rna_qc_passed, genotype_qc_passed, study, protocol, sex))
samples = dplyr::filter(meta_imported, rna_qc_passed, genotype_qc_passed) %>%
  dplyr::mutate(type = ifelse(protocol %in% c("total", "poly(A)"), "RNA-seq", "microarray")) %>%
  dplyr::mutate(is_female = ifelse(sex == "female", 1, 0))

#Summarize metadata
study_summaries = dplyr::group_by(samples, study) %>% 
  dplyr::summarise(individual_count = length(unique(genotype_id)), 
                   sample_count = length(unique(sample_id)), 
                   type = type[1], 
                   cell_types = paste(unique(cell_type), collapse = ", "), 
                   conditions = paste(unique(condition), collapse = ", "),
                   female_fraction = sum(is_female)/length(is_female))
write.table(study_summaries, "results/study_metadata_summary.tsv", sep = "\t", quote = F, row.names = F)



