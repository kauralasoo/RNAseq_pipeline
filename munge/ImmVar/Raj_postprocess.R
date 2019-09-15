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

#Make SE and filter
se = makeSummarizedExperimentFromCountMatrix(mat_final, gene_metadata, sample_metadata, assay_name = "exprs")
#Specify valid chromsomes and valid gene types
valid_gene_types = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                     "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                     "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                     "antisense","sense_intronic","sense_overlapping")
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9")

filtered_se = filterSummarizedExperiment(se, filter_rna_qc = TRUE, filter_genotype_qc = TRUE, 
                                         valid_gene_types = valid_gene_types, valid_chromosomes = valid_chromosomes)

pca_res = eQTLUtils::transformSE_PCA(filtered_se, assay_name = "exprs")

filtered_mat = mat_final[,c("phenotype_id",colnames(filtered_se))]
filtered_mat = dplyr::filter(filtered_mat, phenotype_id %in% rownames(filtered_se))

#Export matrix
write.table(filtered_mat, "results/Raj_2014.hugene_10_ST.filtered.tsv", sep = "\t", quote = FALSE, row.names = FALSE)





