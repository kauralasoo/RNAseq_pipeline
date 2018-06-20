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


#GEUVADIS
read_counts = readr::read_tsv("processed/GEUVADIS/matrices/gene_expression_featureCounts.txt.gz")
sample_metadata = readr::read_delim("metadata/GEUVADIS/GEUVADIS_compiled_metadata.txt", delim = "\t") %>%
  dplyr::mutate(cell_type = "LCL", condition = "naive", read_length = "75bp", stranded = FALSE, paired = TRUE, protocol = "poly(A)") %>%
  dplyr::select(-fq1, -fq2)
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")

#Make GEUVADIS SE object
featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)
processed_se = filterSE_gene_types(featureCounts_se) %>% normaliseSE_cqn()

pca_res = transformSE_PCA(processed_se, n_pcs = 10, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = population_code)) + geom_point()



# BLUEPRINT ---------------------------------------------------------------
read_counts = readr::read_tsv("processed/BLUEPRINT/gene_expression_featureCounts.txt.gz") %>%
  dplyr::filter(!(gene_id %like% "PAR_Y"))
read_counts$gene_id = (dplyr::select(read_counts, gene_id) %>% tidyr::separate(gene_id, c("gene_id", "suffix"), sep = "\\."))$gene_id
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")
sample_metadata = readr::read_csv("metadata/BLUEPRINT/BLUEPRINT_extracted_metadata.csv") %>%
  dplyr::transmute(sample_id = X1, original_cell_type = CELL_TYPE, donor_age = DONOR_AGE, sex = gender, protocol = MOLECULE, genotype_id = donor_id) %>%
  dplyr::mutate(cell_type = case_when(
    original_cell_type == "CD4-positive, alpha-beta T cell" ~ "T-cell",
    original_cell_type == "mature neutrophil" ~ "neutrophil",
    original_cell_type == "CD14-positive, CD16-negative classical monocyte" ~ "monocyte"
  )) %>%
  dplyr::mutate(condition = "naive", read_length = "100bp", stranded = TRUE, protocol) %>%
  dplyr::mutate(protocol = ifelse(protocol == "total RNA", "total", "poly(A)"))
is_paired = readr::read_delim("metadata/BLUEPRINT/PE_snakemake_string.txt", delim = ": ", col_names = c("sample_id", "blah")) %>%
  dplyr::select(sample_id) %>%
  dplyr::mutate(paired = TRUE)
sample_metadata = dplyr::left_join(sample_metadata, is_paired, by = "sample_id") %>% dplyr::mutate(paired = ifelse(is.na(paired), FALSE, paired))

featureCounts_se = makeSummarizedExperiemnt(read_counts, transcript_meta, sample_metadata)
processed_se = filterSE_gene_types(featureCounts_se, valid_gene_types = "protein_coding") %>% normaliseSE_cqn()

pca_res = transformSE_PCA(processed_se, n_pcs = 10, center = TRUE, scale. = TRUE)
ggplot(pca_res$pca_matrix, aes(x = PC1, y = PC2, color = protocol)) + geom_point()

