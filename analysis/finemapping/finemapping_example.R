library("SNPRelate")
library("GDSArray")
library("devtools")
library("ggplot2")
library("susieR")
load_all("../eQTLUtils/")

importQtlmapCovariates <- function(covariates_path){
  pc_matrix = read.table(covariates_path, check.names = F, header = T, stringsAsFactors = F)
  pc_transpose = t(pc_matrix[,-1])
  colnames(pc_transpose) = pc_matrix$SampleID
  pc_df = dplyr::mutate(as.data.frame(pc_transpose), genotype_id = rownames(pc_transpose)) %>%
    dplyr::as_tibble() %>% 
    dplyr::select(genotype_id, dplyr::everything())
  
  #Make PCA matrix
  pc_matrix = as.matrix(dplyr::select(pc_df,-genotype_id))
  rownames(pc_matrix) = pc_df$genotype_id
  return(pc_matrix)
}

#Set parameters
cis_distance = 500000
gds_file = "results/finemapping/genotypes/Alasoo_2018_GRCh38.filtered.gds"
pca_covariates_file = "results/finemapping/Alasoo_2018/PCA/Alasoo_2018_ge_macrophage_naive/Alasoo_2018_ge_macrophage_naive.covariates.txt"

#Make a GDS matrix
SNPRelate::snpgdsVCF2GDS("results/finemapping/genotypes/Alasoo_2018_GRCh38.filtered.vcf.gz", "results/finemapping/genotypes/Alasoo_2018_GRCh38.filtered.gds", method="biallelic.only")

#Make a SummarizedExperiment of the expression data
count_matrix = readr::read_tsv("results/finemapping/Alasoo_2018.gene_counts_cqn_norm.tsv")
sample_metadata = utils::read.csv("../SampleArcheology/studies/cleaned/Alasoo_2018.tsv", sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim("~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz", delim = "\t", col_types = "ccccciiicciidi")
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = count_matrix, 
                                                         row_data = phenotype_meta, 
                                                         col_data = sample_metadata, 
                                                         quant_method = "gene_counts")

#Import variant info
variant_info = importVariantInformationFromGDS(gds_file)

#Extract phenotype from SE
gene_vector = extractPhentypeFromSE("ENSG00000170458", se, "counts") %>%
  dplyr::mutate(phenotype_value_std = (phenotype_value - mean(phenotype_value))/sd(phenotype_value))
gene_meta = dplyr::filter(SummarizedExperiment::rowData(se) %>% as.data.frame(), phenotype_id == "ENSG00000170458")

#Import PCs
pc_matrix = importQtlmapCovariates(pca_covariates_file)
pc_matrix = pc_matrix[naive_df$genotype_id,]

#Import genotype matrix
genotype_matrix = extractGenotypeMatrixFromGDS(
  chr = gene_meta$chromosome, 
  start = gene_meta$phenotype_pos - cis_distance, 
  end = gene_meta$phenotype_pos + cis_distance, 
  variant_information = variant_info, 
  gdsfile = gds_file)

#Residualise gene expression
model_fit = lm.fit(pc_matrix, naive_df$phenotype_value_std)
naive_residuals = dplyr::mutate(naive_df, phenotype_residual = model_fit$residuals) %>%
  dplyr::mutate(phenotype_residual_std = (phenotype_residual - mean(phenotype_residual))/sd(phenotype_residual))

#Fit finemapping model
expression_vector = naive_residuals$phenotype_residual_std
names(expression_vector) = naive_residuals$genotype_id
gt_matrix = genotype_matrix[,names(expression_vector)]
gt_std = t(gt_matrix - apply(gt_matrix, 1, mean))

fitted <- susieR::susie(gt_std, expression_vector,
                        L = 10,
                        estimate_residual_variance = TRUE, 
                        estimate_prior_variance = FALSE,
                        scaled_prior_variance = 0.1,
                        verbose = TRUE,
                        compute_univariate_zscore = TRUE)

#Add variant id vector to the model object
fitted$variant_ids = colnames(gt_std)

susieR::susie_plot(fitted,y="z")
susieR::susie_plot(fitted,y="PIP")


