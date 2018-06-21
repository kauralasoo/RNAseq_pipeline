library("readr")
library("dplyr")
library("devtools")
load_all("../seqUtils/")
library("tidyr")
library("ggplot2")
library("cqn")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")


BLUEPRINT = readRDS("results/SummarizedExperiments/BLUEPRINT.rds")
Quach = readRDS("results/SummarizedExperiments/Quach_2016_Monocytes.rds")
Nedelec = readRDS("results/SummarizedExperiments/Nedelec_2016_Macrophages.rds")
Alasoo = readRDS("results/SummarizedExperiments/Alasoo_2018_Macrophages.rds")
GEUVADIS = readRDS("results/SummarizedExperiments/GEUVADIS.rds")
TwinsUK = readRDS("results/SummarizedExperiments/TwinsUK.rds")


#Merge studies
merged_data = mergeCountsSEs(BLUEPRINT, Quach) %>%
  mergeCountsSEs(Nedelec) %>%
  mergeCountsSEs(Alasoo) %>%
  mergeCountsSEs(GEUVADIS) %>%
  mergeCountsSEs(TwinsUK)


#Keep only naive samples
merged_data = merged_data[,merged_data$condition == "naive"]


processed_se = filterSE_gene_types(merged_data, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = cell_type, shape = study)) + geom_point()

#Perform MDS
matrix = log(assays(processed_se)$tpms+0.1,2)
dist = cor(matrix, method = "pearson")
fit <- MASS::isoMDS(1-dist, k=2)

mds_matrix = as.data.frame(fit$points) %>%
  as_tibble() %>%
  dplyr::mutate(sample_id = rownames(fit$points)) %>%
  dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id")
mds_plot = ggplot(mds_matrix, aes(x = V1, y = V2, color = cell_type, shape = study)) + 
  geom_point() +
  xlab("MDS Coordinate 1") + 
  ylab("MDS Coordinate 2")
ggsave("results/figures/MDS_cell_types_plot.pdf", plot = mds_plot, width = 5, height = 4)


#Explore stimulated conditions
merged_data = merged_data[,merged_data$cell_type %in% c("macrophage")]

processed_se = filterSE_gene_types(merged_data, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition, shape = study)) + geom_point()

