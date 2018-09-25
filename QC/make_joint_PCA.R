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
library("Rtsne")


se_list = list(BLUEPRINT = readRDS("results/SummarizedExperiments/BLUEPRINT.rds"),
               Quach = readRDS("results/SummarizedExperiments/Quach_2016_Monocytes.rds"),
               Nedelec = readRDS("results/SummarizedExperiments/Nedelec_2016_Macrophages.rds"),
               Alasoo = readRDS("results/SummarizedExperiments/Alasoo_2018_Macrophages.rds"),
               GEUVADIS = readRDS("results/SummarizedExperiments/GEUVADIS.rds"),
               TwinsUK = readRDS("results/SummarizedExperiments/TwinsUK_old.rds"),
               Fairfax = readRDS("results/SummarizedExperiments/Fairfax_2018.rds"),
               GENCORD = readRDS("results/SummarizedExperiments/GENCORD.rds"))
#se_list$TwinsUK = se_list$TwinsUK[,se_list$TwinsUK$cell_type == "LCL"]

#Merge studies
merged_data = purrr::reduce(se_list, mergeCountsSEs)

#Exclude low quality samples
TwinsUK_exclude = c("TWPID10593_B", "TWPID11605_B", "TWPID2140_F","TWPID9396_F", "TWPID12889_S", "TWPID8405_S")
GENCORD_exclude = c("UCB1193", "UCF1183", "UCF1119")
merged_data = merged_data[,!(merged_data$sample_id %in% c(TwinsUK_exclude, GENCORD_exclude))]

#Count donors and samples
coldata = colData(merged_data) %>% tbl_df2()

n_samples = dplyr::select(coldata, study, sample_id) %>% 
  dplyr::distinct() %>% dplyr::group_by(study) %>% 
  dplyr::summarise(n_samples = length(sample_id))

n_donors = dplyr::select(coldata, study, genotype_id) %>% 
  dplyr::distinct() %>% dplyr::group_by(study) %>% 
  dplyr::summarise(n_donors = length(genotype_id))

counts = dplyr::left_join(n_samples, n_donors, by = "study")
write.table(counts, "results/counts.txt", sep = "\t", quote = FALSE, row.names = FALSE)



#Keep only naive samples
merged_data = merged_data[,merged_data$condition == "naive"]

#TPM normalisation
processed_se = filterSE_gene_types(merged_data, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()

#Perform MDS
matrix = log(assays(processed_se)$tpms+0.1,2)
dist = cor(matrix, method = "pearson")
fit <- MASS::isoMDS(1-dist, k=2)

mds_matrix = as.data.frame(fit$points) %>%
  as_tibble() %>%
  dplyr::mutate(sample_id = rownames(fit$points)) %>%
  dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id") %>%
  dplyr::mutate(study_label = as.character(as.numeric(study)))
mds_plot = ggplot(mds_matrix, aes(x = V1, y = V2, color = cell_type, label = study_label)) + 
  geom_text() +
  xlab("MDS Coordinate 1") + 
  ylab("MDS Coordinate 2")
ggsave("results/figures/MDS_cell_types_plot.pdf", plot = mds_plot, width = 7, height = 6)


#Make PCA
filtered_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(filtered_se, assay_name = "tpms", n_pcs = NULL, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = study)) + geom_point()



#Explore stimulated conditions
merged_data = merged_data[,merged_data$cell_type %in% c("macrophage")]

processed_se = filterSE_gene_types(merged_data, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = condition, shape = study)) + geom_point()




#Perform UMAP with the uwot package
matrix = log(assays(processed_se)$tpms+0.1,2)

#Use eucledian distance
umap_mat = t(matrix)
umap_res <- uwot::umap(umap_mat, n_neighbors = 50, alpha = 0.5, init = "random")
umap_df = data_frame(umap_1 = umap_res[,1], umap_2 = umap_res[,2], sample_id = rownames(umap_mat)) %>%
  dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id")

#Use pearson correlation
umap_dist = 1-dist
umap_res <- uwot::umap(umap_dist, n_neighbors = 100, alpha = 0.5, init = "pca")
umap_df = data_frame(umap_1 = umap_res[,1], umap_2 = umap_res[,2], sample_id = rownames(umap_dist)) %>%
  dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id")


ggplot(umap_df, aes(x = umap_1, y = umap_2, color = cell_type)) + 
  geom_point()



pca_1_10 = pca_res$pca_matrix[,1:20]
tsne_res = Rtsne(pca_1_10, preplexity = 100, pca = FALSE)

tsne_df = data_frame(dim_1 = tsne_res$Y[,1], dim_2 = tsne_res$Y[,2], sample_id = pca_res$pca_matrix$sample_id) %>%
  dplyr::left_join(as.data.frame(colData(processed_se)), by = "sample_id")

ggplot(tsne_df, aes(x = dim_1, y = dim_2, color = cell_type)) + 
  geom_point()



#Perform PCA for TwinsUK
processed_data = processed_se[,processed_se$study == "TwinsUK"]

#Make PCA
pca_res = transformSE_PCA(processed_data, assay_name = "tpms", n_pcs = NULL, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC3, y = PC4, color = cell_type, label = sample_id)) + geom_point()


#Perform PCA on the GENCORD dataset
gencord = merged_data[,merged_data$study == "GENCORD"]

processed_se = filterSE_gene_types(gencord, valid_gene_types = "protein_coding") %>% normaliseSE_tpm()
processed_se = processed_se[apply(assays(processed_se)$tpms, 1, median) > 1,]
pca_res = transformSE_PCA(processed_se, assay_name = "tpms", n_pcs = 10, log_transform = TRUE, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = cell_type, shape = study, label = sample_id)) + geom_text()

mat = assays(processed_se)$tpms
dist = cor(mat, method = "spearman")

pdf("results/figures/GENCORD_heatmap.pdf", width = 35, height = 35)
gplots::heatmap.2(dist, trace = "none")
dev.off()

TwinsUK_old = readRDS("results/SummarizedExperiments/TwinsUK_old.rds")
TwinsUK_new = readRDS("results/SummarizedExperiments/TwinsUK.rds")

