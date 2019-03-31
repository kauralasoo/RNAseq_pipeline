library("dplyr")
library("devtools")
load_all("../eQTLUtils/")

#Import data
gene_meta = read.table("metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz", header = T, stringsAsFactors = F) %>% dplyr::as_tibble()
sample_metadata = read.table("../SampleArcheology/studies/cleaned/GENCORD.tsv", header = T, stringsAsFactors = F) %>% dplyr::as_tibble()
featureCounts_matrix = read.table("results/expression_matrices/featureCounts/GENCORD.tsv.gz", sep = "\t")

#Make a summarized experiment
featureCounts_se = eQTLUtils::makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "counts")

#Keep high quality samples
filtered_se = eQTLUtils::filterSummarizedExperiment(featureCounts_se, filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
filtered_meta = SummarizedExperiment::colData(filtered_se) %>% as.data.frame() %>% dplyr::as_tibble()


#### Calculate median tpm per qtl_group ####
median_tpm_df = eQTLUtils::estimateMedianTPM(filtered_se, subset_by = "qtl_group", assay_name = "counts")
gzfile = gzfile("results/QC_report/GENCORD/median_tpm_by_qtl_group.txt.gz", "w")
write.table(median_tpm_df, gzfile, sep = "\t", row.names = F, quote = F)
close(gzfile)


## QC steps ###

#### 1. Find best matches with mbv ####
mbv_best_matches = eQTLUtils::mbvImportData("processed/GENCORD/mbv/") %>%
  purrr::map_df(., eQTLUtils::mbvFindBestMatch, .id = "sample_id")
matched_ids = dplyr::select(filtered_meta, sample_id, genotype_id) %>% 
  dplyr::left_join(mbv_best_matches, by = "sample_id") %>%
  dplyr::mutate(is_match = ifelse(genotype_id == mbv_genotype_id, TRUE, FALSE))

#Find and report all mismatches
mismatched_ids = dplyr::filter(matched_ids, !is_match)
if(nrow(mismatched_ids) > 0){
  message("Following genotype mismatches detected:")
  mismatched_ids
}

#### 2. Perform gene expression PCA and export plots ####
#TODO: Colour points according to qtl_group

#### 3. Perform genotype PCA and export plot ####

#### 4. Look for related individuals ####

#### 5. Check the exporession of XIST vs Y-genes and plot those against annotated sex ####
#TODO: Can we automatically flag outliers here?



