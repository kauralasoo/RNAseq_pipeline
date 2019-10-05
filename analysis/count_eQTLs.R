library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
load_all("../eQTLUtils/")

#Import study metadata
meta_list = list(BLUEPRINT_PE = "../SampleArcheology/studies/cleaned/BLUEPRINT_PE.tsv",
                 BLUEPRINT_SE = "../SampleArcheology/studies/cleaned/BLUEPRINT_SE.tsv",
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
  dplyr::mutate(type = ifelse(protocol %in% c("total", "poly(A)"), "RNA-seq", "microarray"))

sample_sizes = dplyr::group_by(samples, study, qtl_group, type) %>% 
  dplyr::summarise(sample_size = n())

#Import RNA-seq gene expression results
gene_meta = read.table("~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz", stringsAsFactors = F, header = TRUE) %>%
  dplyr::as_tibble() %>%
  dplyr::select(phenotype_id, group_id, gene_id, gene_name)
sizes = dplyr::filter(sample_sizes, type == "RNA-seq") %>% 
  dplyr::mutate(file_path = file.path("results/lead_variants/", study, paste0(study, "_ge_", qtl_group, ".lead_variants.txt")))

#Import QTL results
file_list = setNames(as.list(sizes$file_path), sizes$file_path)
qtl_results = purrr::map(file_list, ~eQTLUtils::importQTLtoolsNominalTable(.) %>% dplyr::left_join(gene_meta, by = "phenotype_id"))
qtl_counts = purrr::map_df(qtl_results, ~dplyr::mutate(., n_genes = n()) %>% 
                             dplyr::filter(p_fdr < 0.05) %>% 
                             dplyr::mutate(n_qtls = n()) %>%
                             dplyr::select(n_genes, n_qtls) %>%
                             dplyr::distinct(), .id = "file_path")

#Merge
eqtl_results = dplyr::left_join(sizes, qtl_counts)


#### Import array results ####
gene_meta = read.table("~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz", stringsAsFactors = F, header = TRUE) %>%
  dplyr::as_tibble() %>%
  dplyr::select(phenotype_id, group_id, gene_id, gene_name)
sizes = dplyr::filter(sample_sizes, type == "microarray") %>% 
  dplyr::mutate(file_path = file.path("results/lead_variants/", study, paste0(qtl_group, ".lead_variants.txt"))) %>%
  dplyr::filter(!(qtl_group %in% c("ileum", "transverse_colon", "rectum")))

#Import QTL results
file_list = setNames(as.list(sizes$file_path), sizes$file_path)
qtl_results = purrr::map(file_list, ~eQTLUtils::importQTLtoolsNominalTable(.) %>% dplyr::left_join(gene_meta, by = "phenotype_id"))
qtl_counts = purrr::map_df(qtl_results, ~dplyr::group_by(., gene_id) %>% 
                             dplyr::arrange(gene_id, p_fdr) %>% 
                             dplyr::filter(row_number() == 1) %>%
                             dplyr::ungroup() %>%
                             dplyr::mutate(n_genes = n()) %>%
                             dplyr::filter(p_fdr < 0.05) %>% 
                             dplyr::mutate(n_qtls = n()) %>%
                             dplyr::select(n_genes, n_qtls) %>%
                             dplyr::distinct(), .id = "file_path")

#Merge
array_results = dplyr::left_join(sizes, qtl_counts)
full_results = dplyr::bind_rows(eqtl_results, array_results) %>%
  dplyr::filter(!(qtl_group %in% c("CD8_T-cell_anti-CD3-CD28", "CD4_T-cell_anti-CD3-CD28")))

eGene_scatter = ggplot(full_results, aes(x = sample_size, y = n_qtls, color = study, label = qtl_group)) + 
  geom_point() +
  scale_y_continuous(limits = c(0, 7500)) +
  scale_x_continuous(limits = c(0, 600)) +
  ylab("Number of eGenes") +
  xlab("Sample size") +
  theme_light() + 
  theme(legend.position = "none")
ggsave("results/figures/eGenes_samples_scatter.pdf", plot = eGene_scatter, width= 4, height = 4)



