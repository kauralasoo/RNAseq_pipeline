library("devtools")
library("dplyr")
load_all("../eQTLUtils/")


#Fairfax 2014
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/Fairfax_2014.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz", sep = "\t")
fairfax_2014_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")
fairfax_se_filtered = fairfax_2014_se[,fairfax_2014_se$rna_qc_passed & fairfax_2014_se$genotype_qc_passed]

#Normalize
fairfax_norm = array_normaliseSE(fairfax_se_filtered, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)

#Export expression data to disk
studySEtoQTLTools(fairfax_norm, assay_name = "norm_exprs", "processed/Fairfax_2014/qtltools/input/array/")

#Perform PCA
pca = eQTLUtils::transformSE_PCA(fairfax_norm, assay_name = "norm_exprs")
ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()

