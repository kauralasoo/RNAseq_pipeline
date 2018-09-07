library("readr")
library("devtools")
load_all("../seqUtils/")
library("tidyr")
library(ggplot2)

read_counts = readr::read_tsv("processed/Fairfax/matrices/gene_expression_featureCounts.txt.gz")
gene_lengths = dplyr::select(read_counts, gene_id, length)

sample_names = colnames(read_counts)[3:ncol(read_counts)]

sample_meta = dplyr::data_frame(sample_id = sample_names) %>%
  tidyr::separate(sample_id, c("donor", "condition"), sep = "_", extra = "merge", remove = FALSE) %>%
  dplyr::mutate(read_length = "100bp", stranded = TRUE, paired = TRUE, cell_type = "monocytes")

#Make count matrix
count_matrix = as.matrix(dplyr::select(read_counts, -gene_id, -length))
rownames(count_matrix) = read_counts$gene_id

tpm_matrix = calculateTPM(count_matrix, gene_lengths)
mean_tpm = calculateMean(tpm_matrix, as.data.frame(sample_meta), factor = "condition")
expressed_tpm = tpm_matrix[apply(mean_tpm, 1, min) > 1,]

#Perform PCA
pca_mat = performPCA(log(expressed_tpm + 1,2), sample_meta)
pca_matrix = dplyr::filter(pca_mat$pca_matrix, condition %in% c("UT", "IFN", "LPS24", "MET1m"))
ggplot(pca_matrix, aes(x = PC1, y = PC2, color = condition)) + geom_point()

