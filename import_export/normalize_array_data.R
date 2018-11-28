library("devtools")
library("dplyr")
library("SummarizedExperiment")
load_all("../eQTLUtils/")

#### Fairfax_2014 ####
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/Fairfax_2014.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz", sep = "\t")
fairfax_2014_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

fairfax_2014_norm = filterSummarizedExperiment(fairfax_2014_se, filter_rna_qc = TRUE, filter_genotype_qc = TRUE) %>%
  array_normaliseSE(norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, 
                    adjust_batch = TRUE, filter_quality = TRUE)
saveRDS(fairfax_2014_norm, "results/SummarizedExperiments/array_norm/Fairfax_2014.rds")

#### Fairfax 2012 ####
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/Fairfax_2012.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/Fairfax_2012.tsv.gz", sep = "\t")
fairfax_2012_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

#Filter SE to keep correct chromosomes and QC-passed samples
fairfax_2012_norm = filterSummarizedExperiment(fairfax_2012_se, filter_rna_qc = TRUE, filter_genotype_qc = TRUE) %>%
  array_normaliseSE(norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)
saveRDS(fairfax_2012_norm, "results/SummarizedExperiments/array_norm/Fairfax_2012.rds")


#### CEDAR ####
gene_meta = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz")
sample_metadata = readr::read_tsv("metadata/cleaned/CEDAR.tsv")
expression_matrix = read.table("results/expression_matrices/HumanHT-12_V4/CEDAR.tsv.gz", sep = "\t", check.names = FALSE)
cedar_se = makeSummarizedExperiment(expression_matrix, gene_meta, sample_metadata, assay_name = "exprs")

#Filter SE to keep correct chromosomes and QC-passed samples
cedar_norm = filterSummarizedExperiment(cedar_se, filter_rna_qc = TRUE, filter_genotype_qc = TRUE) %>%
  array_normaliseSE(norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = TRUE)
saveRDS(cedar_norm, "results/SummarizedExperiments/array_norm/CEDAR.rds")

