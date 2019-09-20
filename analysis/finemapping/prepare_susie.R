library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
load_all("../eQTLUtils/")

makeSusieInputTable <- function(qtl_groups_df){
  
  #quant methods
  quants = c("ge","exon","tx", "txrev")
  quant_method = c("gene_counts", "exon_counts", "transcript_usage", "txrevise")
  matrix_suffix = c("gene_counts_cqn_norm.tsv","exon_counts_cqn_norm.tsv","transcript_usage_qnorm.tsv", "txrevise_qnorm.tsv")
  meta_file_name = c("gene_counts_Ensembl_96_phenotype_metadata.tsv.gz","exon_counts_Ensembl_96_phenotype_metadata.tsv.gz","transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz","txrevise_Ensembl_96_phenotype_metadata.tsv.gz")
  quant = dplyr::tibble(quant_method = quant_method, quant = quants, matrix_suffix = matrix_suffix, meta_file_name = meta_file_name, c = 1)
  
  #Input paths
  sumstats_home = "/gpfs/hpc/home/a72094/datasets/summary_stats/eQTLCatalogue/v0.1/"
  matrix_home = "/gpfs/hpc/home/a72094/datasets/processed/expression_matrices/normalised_studies/"
  sample_meta_home = "/gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/studies/cleaned/"
  phenotype_meta_home = "/gpfs/hpc/home/a72094/annotations/eQTLCatalogue/v0.1/phenotype_metadata/"
  
  #Add uqant method to df
  qtl_df = dplyr::mutate(qtl_groups_df, c = 1) %>%
    dplyr::left_join(quant, by = "c") %>%
    dplyr::select(-c) %>%
    dplyr::mutate(set = paste(study, quant, qtl_group, sep = "_")) %>%
    dplyr::mutate(covariates = file.path(sumstats_home, "pipeline_output", study, "PCA", set, paste0(set, ".covariates.txt"))) %>%
    dplyr::mutate(phenotype_list = file.path(sumstats_home, "summary_stats", study,paste0(set,".nominal.sorted.txt.gz"))) %>%
    dplyr::mutate(expression_matrix = file.path(matrix_home, study, paste0(study, ".", matrix_suffix))) %>%
    dplyr::mutate(sample_meta = file.path(sample_meta_home, paste0(study, ".tsv"))) %>%
    dplyr::mutate(phenotype_meta = file.path(phenotype_meta_home, meta_file_name)) %>%
    dplyr::select(study, qtl_group, quant_method, expression_matrix, phenotype_meta, sample_meta, vcf, phenotype_list, covariates)
  return(qtl_df)
}

#Import metadata files
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
                 Ye_2018 = "../SampleArcheology/studies/cleaned/Ye_2018.tsv",
                 Raj_2014 = "../SampleArcheology/studies/cleaned/Raj_2014.tsv",
                 Schmiedel_2018 = "../SampleArcheology/studies/cleaned/Schmiedel_2018.tsv",
                 Lepik_2017 = "../SampleArcheology/studies/cleaned/Lepik_2017.tsv",
                 ROSMAP = "../SampleArcheology/studies/cleaned/ROSMAP.tsv",
                 BrainSeq = "../SampleArcheology/studies/cleaned/BrainSeq.tsv")
meta_imported = purrr::map_df(meta_list, ~read.table(., sep = "\t", stringsAsFactors = F, header =T) %>% 
                                dplyr::filter(genotype_qc_passed, rna_qc_passed) %>%
                                dplyr::as_tibble() %>%
                                dplyr::select(sample_id, genotype_id, cell_type, condition, qtl_group, rna_qc_passed, genotype_qc_passed, study, protocol, sex))

#Import VCF file paths
vcf_paths = read.table("~/projects/SampleArcheology/studies/study_genotype_paths.tsv", stringsAsFactors = F, header = TRUE) %>%
  dplyr::as_tibble()

#Extract studies and qtl groups
qtl_groups_df = dplyr::filter(meta_imported, study %in% 
                                c("BLUEPRINT_PE","BLUEPRINT_SE","Alasoo_2018", "Quach_2016", "Nedelec_2016", 
                                  "Schmiedel_2018","GEUVADIS", "GENCORD", "TwinsUK", "Lepik_2017")) %>% 
  dplyr::transmute(study, qtl_group) %>% 
  dplyr::distinct() %>% 
  dplyr::left_join(vcf_paths, by = "study")

susie_file = makeSusieInputTable(qtl_groups_df) %>% 
  dplyr::filter(susie_file, quant_method != "exon_counts")
write.table(susie_file, "../SampleArcheology/finemapping/Alasoo_2018.tsv", sep = "\t", quote = F, row.names = F)
