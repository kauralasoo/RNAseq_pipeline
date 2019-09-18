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

#Array datasets
array_study_list = dplyr::filter(meta_imported, study %in% c("CEDAR", "Fairfax_2014", "Fairfax_2012", "Kasela_2017", "Naranbhai_2015")) %>% 
  dplyr::transmute(study, qtl_group, quant_method = protocol) %>% 
  dplyr::distinct() %>%
  dplyr::mutate(qtl_leads = file.path("/gpfs/hpc/home/a72094/datasets/summary_stats/eQTLCatalogue/v0.1/summary_stats/CEDAR/",paste0(qtl_group, ".permuted.txt.gz"))) %>%
  dplyr::mutate(qtl_stats = file.path("/gpfs/hpc/home/a72094/datasets/summary_stats/eQTLCatalogue/v0.1/summary_stats/CEDAR/",paste0(qtl_group, ".nominal.sorted.txt.gz"))) %>%
  dplyr::mutate(qtl_varinfo = file.path("/gpfs/hpc/home/a72094/datasets/summary_stats/eQTLCatalogue/v0.1/summary_stats/CEDAR/",paste0(qtl_group, ".variant_information.txt.gz")))

write.table(array_study_list, "../colocWrapper/testdata/study_file.tsv", sep = "\t", row.names = F, quote = F)

