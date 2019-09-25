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
                 #Ye_2018 = "../SampleArcheology/studies/cleaned/Ye_2018.tsv",
                 #Raj_2014 = "../SampleArcheology/studies/cleaned/Raj_2014.tsv",
                 Schmiedel_2018 = "../SampleArcheology/studies/cleaned/Schmiedel_2018.tsv",
                 ROSMAP = "../SampleArcheology/studies/cleaned/ROSMAP.tsv",
                 BrainSeq = "../SampleArcheology/studies/cleaned/BrainSeq.tsv",
                 Lepik_2017 = "../SampleArcheology/studies/cleaned/Lepik_2017.tsv")
meta_imported = purrr::map_df(meta_list, ~read.table(., sep = "\t", stringsAsFactors = F, header =T) %>% 
                             dplyr::as_tibble() %>%
                             dplyr::transmute(sample_id, genotype_id, cell_type, condition, qtl_group, rna_qc_passed, genotype_qc_passed, study, protocol, sex, timepoint = as.character(timepoint)))
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

#Import ontology mapping
mappings = read.table("../eQTL-Catalogue-resources/ontology_mappings/tissue_onotology_mapping.tsv", header =T, sep = "\t", stringsAsFactors = F)
friendly_names = read.table("../eQTL-Catalogue-resources/ontology_mappings/friendly_names.tsv", sep = "\t", header = TRUE, stringsAsFactors = F)
levels = c("monocyte", "CD16+ monocyte", "macrophage", "dendritic cell", "neutrophil", "platelet","CD4+ T cell", "CD8+ T cell", 
           "T cell", "Tfh cell","Th17 cell","Th1 cell","Th2 cell","Treg naive","Treg memory","NK cell", "B cell", "LCL","blood","DLPFC",
           "fibroblast","adipose","skin","pancreatic islet", "iPSC","sensory neuron","transverse colon", "ileum", "rectum")
friendly_names$ontology_tissue = factor(friendly_names$ontology_tissue, levels = levels)

#Count unique tissues
tissue_counts = dplyr::left_join(samples, mappings, by = c("study", "qtl_group","cell_type")) %>%
  dplyr::left_join(friendly_names, by = c("ontology_term", "ontology_label")) %>%
  dplyr::filter(condition %in% c("memory","naive")) %>%
  dplyr::filter(!(qtl_group %in% c("Th1-17_memory"))) %>%
  dplyr::select(ontology_tissue, study, type) %>%
  dplyr::group_by(ontology_tissue, study, type) %>%
  dplyr::summarise(n = n())

#RNA-seq datasets
seq_plot = ggplot(tissue_counts, aes(x= ontology_tissue, y = n, fill = study)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Sample size") + 
  xlab("Cell type or tissue")
ggsave("results/figures/sample_stats/sample_counts.pdf", plot = seq_plot, width = 6, height = 3)

#Count conditions
condition_mappings = read.table("../eQTL-Catalogue-resources/ontology_mappings/condition_ontology_mapping.tsv", header =T, sep = "\t", stringsAsFactors = F)
condition_counts = dplyr::filter(samples, !(condition %in% c("memory","naive"))) %>%
  dplyr::left_join(mappings, by = c("study", "qtl_group","cell_type")) %>%
  dplyr::left_join(friendly_names, by = c("ontology_term", "ontology_label")) %>%
  dplyr::select(study, condition, ontology_tissue, type) %>%
  dplyr::filter(condition != "LPS24") %>%
  dplyr::left_join(condition_mappings, by = c("study","condition")) %>%
  dplyr::group_by(ontology_condition, ontology_tissue) %>%
  dplyr::summarise(n = n())

cond_plot = ggplot(condition_counts, aes(x= ontology_condition, y = n, fill = ontology_tissue)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Sample size") + 
  xlab("Biological context")
ggsave("results/figures/sample_stats/contexts_counts.pdf", plot = cond_plot, width = 5, height = 3)


#Identify cell types and tissues
qtl_groups = dplyr::select(samples, study, cell_type, qtl_group) %>% distinct()
write.table(qtl_groups, "results/study_qtl_groups.tsv", sep = "\t", quote = F, row.names = F)

#Identify cell types and tissues
conditions = dplyr::select(samples, study, condition) %>% distinct()
write.table(conditions, "results/study_conditions.tsv", sep = "\t", quote = F, row.names = F)



