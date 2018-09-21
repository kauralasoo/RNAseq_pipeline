library("readr")
library("dplyr")
library("devtools")
load_all("../eQTLUtils/")
library("tidyr")
library("ggplot2")
library("cqn")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")

#Specify mandatory metadata columns
mandatory_cols = c("sample_id", "genotype_id", "sex", "cell_type", "condition", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")


#GEUVADIS
read_counts = readr::read_tsv("processed/GEUVADIS/matrices/gene_expression_featureCounts.txt") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
mbv_matches = readr::read_tsv("metadata/GEUVADIS/GEUVADIS_mbv_best_match.txt") %>%
  dplyr::select(sample_id, mbv_genotype_id)
sample_metadata = readr::read_delim("metadata/GEUVADIS/GEUVADIS_compiled_metadata.txt", delim = "\t") %>%
  dplyr::mutate(cell_type = "LCL", condition = "naive", read_length = "75bp", stranded = FALSE, paired = TRUE, protocol = "poly(A)", timepoint = 0, rna_qc_passed = TRUE, genotype_qc_passed = TRUE, study = "GEUVADIS") %>%
  dplyr::select(-fq1, -fq2) %>%
  dplyr::select(mandatory_cols, everything()) %>%
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(mbv_genotype_id), FALSE, TRUE)) %>%
  dplyr::select(-mbv_genotype_id)
write.table(sample_metadata, "metadata/cleaned/GEUVADIS.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Make GEUVADIS SE object
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
saveRDS(featureCounts_se, "results/SummarizedExperiments/GEUVADIS.rds")

processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = superpopulation_code)) + geom_point()



# BLUEPRINT ---------------------------------------------------------------
read_counts = readr::read_tsv("processed/BLUEPRINT/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
sample_metadata = readr::read_csv("metadata/BLUEPRINT/BLUEPRINT_extracted_metadata.csv") %>%
  dplyr::transmute(sample_id = X1, original_cell_type = CELL_TYPE, donor_age = DONOR_AGE, sex = gender, protocol = MOLECULE, genotype_id = donor_id) %>%
  dplyr::mutate(cell_type = case_when(
    original_cell_type == "CD4-positive, alpha-beta T cell" ~ "T-cell",
    original_cell_type == "mature neutrophil" ~ "neutrophil",
    original_cell_type == "CD14-positive, CD16-negative classical monocyte" ~ "monocyte"
  )) %>%
  dplyr::mutate(condition = "naive", read_length = "100bp", stranded = TRUE, timepoint = 0, study = "BLUEPRINT") %>%
  dplyr::mutate(protocol = ifelse(protocol == "total RNA", "total", "poly(A)"))
is_paired = readr::read_delim("metadata/BLUEPRINT/PE_snakemake_string.txt", delim = ": ", col_names = c("sample_id", "blah")) %>%
  dplyr::select(sample_id) %>%
  dplyr::mutate(paired = TRUE)
sample_metadata_final = dplyr::left_join(sample_metadata, is_paired, by = "sample_id") %>% 
  dplyr::mutate(paired = ifelse(is.na(paired), FALSE, paired)) %>%
  dplyr::mutate(qtl_mapping = ifelse(protocol == "poly(A)", FALSE, TRUE)) %>%
  dplyr::mutate(qtl_mapping = ifelse(cell_type %in% c("neutrophil","monocyte") & paired, FALSE, qtl_mapping)) %>%
  dplyr::select(mandatory_cols, everything())

#Make SE
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata_final)

#Export data
write.table(sample_metadata_final, "metadata/cleaned/BLUEPRINT.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/BLUEPRINT.rds")

#Perform some QC
processed_se = filterSE_gene_types(featureCounts_se, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = protocol)) + geom_point()



# Nedelec_2016_Macrophages ------------------------------------------------
read_counts = readr::read_tsv("processed/Macrophages_Nedelec_2016/matrices/gene_expression_featureCounts.txt") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
sample_metadata = readr::read_tsv("metadata/Nedelec_2016/Nedelec_2016_compiled_metadata.txt") %>%
  dplyr::transmute(sample_id, genotype_id = donor, sex = "male", cell_type = "macrophage", 
                   condition, timepoint = 5, read_length = "100bp", stranded = FALSE, paired = FALSE, 
                   protocol = "poly(A)", qtl_mapping = TRUE, run_id, geo_sample_id, sequencing_center, flowcell, lane, 
                   self_reported_ethnic, AF_admixture, EU_admixture, ethnic_bin) %>%
  dplyr::mutate(condition = ifelse(condition == "Non-infected", "naive", condition), study = "Nedelec_2016") %>%
  dplyr::select(mandatory_cols, everything())

#Make SE
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
write.table(sample_metadata, "metadata/cleaned/Nedelec_2016_Macrophages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Nedelec_2016_Macrophages.rds")


#QC
processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = sequencing_center)) + geom_point()


# Quach_2016_Monocytes ----------------------------------------------------
read_counts = readr::read_tsv("processed/Monocytes_Quach_2016/matrices/gene_expression_featureCounts.txt") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
outlier_samples = c("LVOMJJWNYB_5")
sample_metadata = readr::read_tsv("metadata/Quach_2016/Quach_2016_compiled_metadata.txt") %>%
  dplyr::transmute(sample_id, genotype_id = donor, sex = "male", cell_type = "monocyte", 
                   condition = condition_name, timepoint = 6, read_length = "100bp", stranded = FALSE, paired = FALSE, 
                   protocol = "poly(A)", qtl_mapping = TRUE, ega_id) %>%
  dplyr::mutate(qtl_mapping = ifelse(sample_id %in% outlier_samples, FALSE, TRUE), study = "Quach_2016") %>%
  dplyr::select(mandatory_cols, everything())
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
write.table(sample_metadata, "metadata/cleaned/Quach_2016_Monocytes.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Quach_2016_Monocytes.rds")


processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + geom_point() + geom_text()




# Alasoo_2018_Macrophages -------------------------------------------------
read_counts = readr::read_tsv("processed/Macrophages_Alasoo_2018/matrices/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
sample_metadata = readr::read_tsv("metadata/Alasoo_2018/RNA_sample_metadata.txt") %>%
  dplyr::mutate(cell_type = "macrophage", condition = case_when(
    condition_name == "naive" ~ "naive",
    condition_name == "IFNg" ~ "IFNg",
    condition_name == "SL1344" ~ "Salmonella",
    condition_name == "IFNg_SL1344" ~ "IFNg+Salmonella")) %>%
  dplyr::mutate(timepoint = case_when(
    condition_name == "naive" ~ 0,
    condition_name == "IFNg" ~ 18,
    condition_name == "SL1344" ~ 5,
    condition_name == "IFNg_SL1344" ~ 5), read_length = "75bp", stranded = TRUE, paired = TRUE, protocol = "poly(A)", qtl_mapping = TRUE, study = "Alasoo_2018") %>%
  dplyr::select(mandatory_cols, everything())
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

write.table(sample_metadata, "metadata/cleaned/Alasoo_2018_Macrophages.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Alasoo_2018_Macrophages.rds")


processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()



# TwinsUK -----------------------------------------------------------------
read_counts = readr::read_tsv("processed/TwinsUK/matrices/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
sample_metadata = dplyr::data_frame(sample_id = colnames(read_counts)[3:ncol(read_counts)]) %>%
  dplyr::mutate(genotype_id = sample_id, sex = as.character(NA), cell_type = "LCL", condition = "naive", timepoint = 0, 
                read_length = "50bp", stranded = FALSE, paired = TRUE, qtl_mapping = TRUE, study = "TwinsUK", protocol = "poly(A)") %>%
  dplyr::select(mandatory_cols, everything())
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

write.table(sample_metadata, "metadata/cleaned/TwinsUK.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/TwinsUK.rds")

processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = cell_type)) + geom_point()


# Fairfax 2018 monocytes --------------------------------------------------
read_counts = readr::read_tsv("processed/Fairfax/matrices/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
mbv_results = readr::read_tsv("../Fairfax_monocytes/data/metadata/Fairfax_mbv_best_match.txt") %>% 
  dplyr::select(sample_id, mbv_genotype_id)
sample_metadata = readr::read_tsv("../Fairfax_monocytes/data/metadata/Fairfax_compiled_metadata.txt") %>%
  dplyr::left_join(mbv_results) %>%
  dplyr::transmute(sample_id, genotype_id = mbv_genotype_id, sex = NA_character_, cell_type = "monocyte", condition_original = condition, 
                   timepoint = 0, read_length, stranded, paired, protocol = "poly(A)", qtl_mapping = TRUE, study = "Fairfax_2018") %>%
  dplyr::mutate(condition = ifelse(condition_original == "UT", "naive", condition_original)) %>%
  dplyr::mutate(timepoint = ifelse(condition == "naive", 0, 24)) %>%
  dplyr::select(mandatory_cols, everything()) %>%
  dplyr::mutate(qtl_mapping = ifelse(condition %in% c("naive", "LPS24", "IFN", "MET1m"), TRUE, FALSE))
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

write.table(sample_metadata, "metadata/cleaned/Fairfax_2018.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Fairfax_2018.rds")

#Keep only relevant conditions
filtered_se = featureCounts_se[,featureCounts_se$qtl_mapping]
processed_se = filterSE_gene_types(filtered_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition, label = sample_id)) + geom_point() + geom_text()




