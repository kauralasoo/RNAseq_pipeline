library(dplyr)
library(SummarizedExperiment)
devtools::load_all("../eQTLUtils/")

#Load the rds object
dataset = readRDS("results/Raj_2014_hugene_eset.rds")

#Import gene metadata
gene_metadata = read.table("~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/Affy_Human_Gene_1_0_ST_Ensembl_96_phenotype_metadata.tsv.gz", stringsAsFactors = F, header = TRUE) %>% dplyr::as_tibble()
sample_metadata = read.table("~/projects/SampleArcheology/studies/cleaned/Raj_2014.tsv", stringsAsFactors = F, header = TRUE, sep = "\t") %>% dplyr::as_tibble()

#Rename rows and columns
mat = exprs(dataset)
rownames(mat) = paste0("AFFY_",rownames(mat))
col_names = dplyr::tibble(sample_string = colnames(mat)) %>% 
  tidyr::separate(sample_string, into = c("sample_id", "rest"), sep = "_", extra = "merge")
colnames(mat) = col_names$sample_id
mat_df = as.data.frame(mat)
mat_df$phenotype_id = rownames(mat_df)
mat_final = dplyr::select(mat_df, phenotype_id, everything()) %>% dplyr::as_tibble()

#Export matrix
write.table(mat_final, "results/Raj_2014.hugene_10_ST.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

se = makeSummarizedExperimentFromCountMatrix(mat_final, gene_metadata, sample_metadata, assay_name = "exprs")
pca_res = eQTLUtils::transformSE_PCA(se,assay_name = "exprs")



