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
mandatory_cols = c("sample_id", "genotype_id", "sex", "cell_type", "condition", "qtl_group", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")
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
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(mbv_genotype_id), FALSE, TRUE)) %>%
  dplyr::select(-mbv_genotype_id) %>%
  dplyr::mutate(qtl_group = cell_type) %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/GEUVADIS.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Make GEUVADIS SE object
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
saveRDS(featureCounts_se, "results/SummarizedExperiments/GEUVADIS.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/GEUVADIS.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)


#Export featureCounts gene metadata
gene_meta = rowData(featureCounts_se)
gz1 = gzfile("metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz","w") 
write.table(gene_meta, gz1, sep = "\t", quote = F, row.names = F)
close(gz1)

#Perform PCA
processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = superpopulation_code)) + geom_point()



# BLUEPRINT ---------------------------------------------------------------
read_counts = readr::read_tsv("processed/BLUEPRINT/matrices/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
mbv_matches = readr::read_tsv("metadata/BLUEPRINT/BLUEPRINT_mbv_best_match.txt") %>%
  dplyr::select(sample_id, mbv_genotype_id)
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
  dplyr::mutate(rna_qc_passed = ifelse(protocol == "poly(A)", FALSE, TRUE)) %>%
  dplyr::mutate(rna_qc_passed = ifelse(cell_type %in% c("neutrophil","monocyte") & paired, FALSE, rna_qc_passed)) %>%
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(mbv_genotype_id), FALSE, TRUE)) %>%
  dplyr::mutate(rna_qc_passed = rna_qc_passed & genotype_qc_passed) %>%
  dplyr::select(-mbv_genotype_id) %>%
  dplyr::mutate(sex = ifelse(genotype_id == "S00PWE", "female", sex)) %>%
  dplyr::mutate(sex = ifelse(genotype_id == "S00PVG", "male", sex)) %>%
  dplyr::mutate(qtl_group = cell_type) %>%
  dplyr::select(mandatory_cols, everything())

#Fix sex
sex_qc = readr::read_tsv("metadata/BLUEPRINT/BLUEPRINT_Sex_QC_table.tsv") %>%
  dplyr::mutate(new_sex = ifelse(log(Y_chrom_mean+1,2) < 1, "female", "male"))

#Make SE
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata_final)

#Export data
write.table(sample_metadata_final, "metadata/cleaned/BLUEPRINT.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/BLUEPRINT.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/BLUEPRINT.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)

#Perform some QC
processed_se = filterSE_gene_types(featureCounts_se, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = protocol)) + geom_point()



# Nedelec_2016 ------------------------------------------------
read_counts = readr::read_tsv("processed/Nedelec_2016/matrices/gene_expression_featureCounts.txt") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
mbv_matches = readr::read_tsv("metadata/Nedelec_2016/Nedelec_2016_mbv_best_match.txt") %>%
  dplyr::select(sample_id, mbv_genotype_id)

sample_metadata = readr::read_tsv("metadata/Nedelec_2016/Nedelec_2016_compiled_metadata.txt") %>%
  dplyr::transmute(sample_id, genotype_id = donor, sex = "male", cell_type = "macrophage", 
                   condition, timepoint = 5, read_length = "100bp", stranded = FALSE, paired = FALSE, 
                   protocol = "poly(A)", rna_qc_passed = TRUE, run_id, geo_sample_id, sequencing_center, flowcell, lane, 
                   self_reported_ethnic, AF_admixture, EU_admixture, ethnic_bin) %>%
  dplyr::mutate(condition = ifelse(condition == "Non-infected", "naive", condition), study = "Nedelec_2016") %>%
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(mbv_genotype_id), FALSE, TRUE)) %>%
  dplyr::select(-mbv_genotype_id) %>%
  dplyr::mutate(qtl_group = paste(cell_type, condition, sep = "_")) %>%
  dplyr::select(mandatory_cols, everything())

#Make SE
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
write.table(sample_metadata, "metadata/cleaned/Nedelec_2016.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Nedelec_2016.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/Nedelec_2016.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)

#QC
processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()


# Quach_2016_Monocytes ----------------------------------------------------
read_counts = readr::read_tsv("processed/Quach_2016/matrices/gene_expression_featureCounts.txt") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
outlier_samples = c("LVOMJJWNYB_5")
mbv_matches = readr::read_tsv("metadata/Quach_2016/Quach_2016_mbv_best_match.txt") %>%
  dplyr::select(sample_id, mbv_genotype_id)
sample_metadata = readr::read_tsv("metadata/Quach_2016/Quach_2016_compiled_metadata.txt") %>%
  dplyr::transmute(sample_id, genotype_id = donor, sex = "male", cell_type = "monocyte", 
                   condition = condition_name, timepoint = 6, read_length = "100bp", stranded = FALSE, paired = FALSE, 
                   protocol = "poly(A)", rna_qc_passed = TRUE, genotype_qc_passed = TRUE, ega_id, study = "Quach_2016") %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %in% outlier_samples, FALSE, TRUE)) %>%
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(genotype_id = mbv_genotype_id) %>%
  dplyr::select(-mbv_genotype_id) %>%
  dplyr::mutate(qtl_group = paste(cell_type, condition, sep = "_")) %>%
  dplyr::select(mandatory_cols, everything())
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
write.table(sample_metadata, "metadata/cleaned/Quach_2016.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Quach_2016.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/Quach_2016.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)

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
  dplyr::mutate(qtl_group = paste(cell_type, condition, sep = "_")) %>%
  dplyr::mutate(rna_qc_passed = TRUE, genotype_qc_passed = TRUE) %>%
  dplyr::select(mandatory_cols, everything())
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

write.table(sample_metadata, "metadata/cleaned/Alasoo_2018.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/Alasoo_2018.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/Alasoo_2018.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)


processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()



# TwinsUK -----------------------------------------------------------------
read_counts = readr::read_tsv("processed/TwinsUK/matrices/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
outlier_samples = c("TWPID10593_B", "TWPID11605_B", "TWPID8405_S","TWPID12889_S", "TWPID2140_F","TWPID9396_F")
mbv_matches = readr::read_tsv("metadata/TwinsUK/TwinsUK_mbv_best_match.txt") %>%
  dplyr::select(sample_id, mbv_genotype_id)
sample_metadata = readr::read_delim("metadata/TwinsUK/TwinsUK_sample_metadata.txt", delim = "\t") %>%
  tidyr::separate(genotype_id, into = c("genotype_id", "suffix"), sep = "_") %>%
  dplyr::select(-suffix) %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %in% outlier_samples, FALSE, TRUE)) %>%
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(id_match = genotype_id == mbv_genotype_id) %>%
  dplyr::mutate(id_match = ifelse(is.na(id_match), FALSE, id_match)) %>%
  dplyr::group_by(mbv_genotype_id, cell_type) %>% 
  dplyr::arrange(mbv_genotype_id, cell_type, -id_match) %>%
  dplyr::mutate(genotype_qc_passed = ifelse(row_number() == 1, TRUE, FALSE)) %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(mbv_genotype_id), FALSE, genotype_qc_passed)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(genotype_id = mbv_genotype_id) %>%
  dplyr::select(-id_match, -mbv_genotype_id, -qtl_mapping) %>%
  dplyr::mutate(qtl_group = cell_type) %>%
  dplyr::select(mandatory_cols, everything())

featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

write.table(sample_metadata, "metadata/cleaned/TwinsUK.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
saveRDS(featureCounts_se, "results/SummarizedExperiments/TwinsUK.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/TwinsUK.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)

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


# GENCORD -----------------------------------------------------------------
read_counts = read.table("processed/GENCORD/matrices/gene_expression_featureCounts.txt", header = T, stringsAsFactors = F) %>%
  dplyr::as_tibble() %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()
mbv_matches = readr::read_tsv("metadata/GENCORD/GENCORD_mbv_best_match.txt") %>%
  dplyr::filter(hom_consistent_frac >= 0.75) %>%
  dplyr::select(sample_id, mbv_genotype_id)
sample_metadata = readr::read_delim("metadata/GENCORD/GENCORD_metadata_final.tsv", delim = "\t") %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %in% c("UCB1193", "UCB087", "UCF1183", "UCF1119"), FALSE, TRUE)) %>%
  dplyr::left_join(mbv_matches, by = "sample_id") %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(mbv_genotype_id), FALSE, TRUE)) %>%
  dplyr::select(-mbv_genotype_id, -qtl_mapping) %>%
  dplyr::mutate(qtl_group = cell_type) %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/GENCORD.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Make GEUVADIS SE object
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Export data
saveRDS(featureCounts_se, "results/SummarizedExperiments/GENCORD.rds")

#Save expression matrix
mat = assays(featureCounts_se)$counts
gz2 = gzfile("results/expression_matrices/featureCounts/GENCORD.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)

#Make PCA
filtered_se = featureCounts_se[,featureCounts_se$rna_qc_passed & featureCounts_se$genotype_qc_passed]
processed_se = filterSE_gene_types(filtered_se) %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = cell_type)) + geom_point()


#HipSci
#Specify mandatory metadata columns
mandatory_cols = c("sample_id", "genotype_id", "sex", "cell_type", "condition", "qtl_group", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")

#GEUVADIS
read_counts = readr::read_tsv("processed/HipSci/matrices/gene_expression_featureCounts.txt") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y")) %>%
  removeGeneVersion()

sample_metadata = readr::read_delim("../SampleArcheology/studies/cleaned/HipSci.tsv", delim = "\t")

#Make SE object
featureCounts_se = makeFeatureCountsSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)

#Filter and normalize for QTL mapping
normalised_se = qtltoolsPrepareSE(featureCounts_se, quant_method = "featureCounts", filter_genotype_qc = FALSE)
cqn_matrix = assays(normalised_se)$cqn
cqn_matrix = round(cqn_matrix, 3)
cqn_df = dplyr::as_tibble(cqn_matrix) %>%
  dplyr::mutate(phenotype_id = rownames(cqn_matrix)) %>%
  dplyr::select(phenotype_id, everything())
write.table(cqn_df, "results/expression_matrices/featureCounts/HipSci.cqn_normalized.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
