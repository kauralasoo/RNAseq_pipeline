library("devtools")
library("dplyr")
library("SummarizedExperiment")
load_all("../eQTLUtils/")


#keep chromosomes
valid_chromosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9")

#Fairfax 2014
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/Fairfax_2014.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz", sep = "\t")
fairfax_2014_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

#Filter SE to keep correct chromosomes and QC-passed samples
fairfax_se_filtered = filterSummarizedExperiment(fairfax_2014_se, valid_chromosomes = valid_chromosomes, filter_rna_qc = TRUE, filter_genotype_qc = TRUE)

#Normalize and regress out batch effects
fairfax_norm = array_normaliseSE(fairfax_se_filtered, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)

#Count the number of variants proximal to each gene and remove genes without variants
var_info = importVariantInformation("processed/Fairfax_2014/qtltools/output/array/final/monocyte_naive.variant_information.txt.gz")
fairfax_norm_filtered = checkCisVariants(fairfax_norm, var_info)

#Export expression data to disk
studySEtoQTLTools(fairfax_norm_filtered, assay_name = "norm_exprs", "processed/Fairfax_2014/qtltools/input/array/")

#Perform PCA
pca = eQTLUtils::transformSE_PCA(fairfax_norm, assay_name = "norm_exprs")
ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()


#### Fairfax 2012 ####
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/Fairfax_2012.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/Fairfax_2012.tsv.gz", sep = "\t")
fairfax_2012_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

#Filter SE to keep correct chromosomes and QC-passed samples
fairfax_se_filtered = filterSummarizedExperiment(fairfax_2012_se, valid_chromosomes = valid_chromosomes, filter_rna_qc = TRUE, filter_genotype_qc = TRUE)

#Normalize and regress out batch effects
fairfax_norm = array_normaliseSE(fairfax_se_filtered, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)

#Count the number of variants proximal to each gene and remove genes without variants
var_info = importVariantInformation("processed/Fairfax_2014/qtltools/output/array/final/monocyte_naive.variant_information.txt.gz")
fairfax_norm_filtered = checkCisVariants(fairfax_norm, var_info)

#Export expression data to disk
studySEtoQTLTools(fairfax_norm_filtered, assay_name = "norm_exprs", "processed/Fairfax_2012/qtltools/input/array/")

#Perform PCA
pca = eQTLUtils::transformSE_PCA(fairfax_norm, assay_name = "norm_exprs")
ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()


#### CEDAR ####
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/CEDAR.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/CEDAR.tsv.gz", sep = "\t", check.names = FALSE)
cedar_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

#Filter SE to keep correct chromosomes and QC-passed samples
se_filtered = filterSummarizedExperiment(cedar_se, valid_chromosomes = valid_chromosomes, filter_rna_qc = TRUE, filter_genotype_qc = TRUE)

#Normalize and regress out batch effects
se_norm = array_normaliseSE(se_filtered, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)

#Count the number of variants proximal to each gene and remove genes without variants
var_info = importVariantInformation("results/var_info/CEDAR_GRCh38.variant_information.txt.gz")
se_norm_filtered = checkCisVariants(se_norm, var_info)

#Export expression data to disk
studySEtoQTLTools(se_norm_filtered, assay_name = "norm_exprs", "processed/CEDAR/qtltools/input/array/")

#Perform PCA
pca = eQTLUtils::transformSE_PCA(se_norm, assay_name = "norm_exprs")
ggplot(pca$pca_matrix, aes(x = PC3, y = PC4, color = cell_type)) + geom_point()

