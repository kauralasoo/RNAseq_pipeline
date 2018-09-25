library("cqn")
library("dplyr")
library("devtools")
library("SummarizedExperiment")
load_all("../seqUtils/")

#Import SE object
se_object = readRDS("results/SummarizedExperiments/GEUVADIS.rds")

#Keep only samples that pass QC
filtered_se = se_object[,se_object$qtl_mapping == TRUE & se_object$superpopulation_code == "EUR"]

#Normalize gene expression
expressed_counts = rowSums(assays(filtered_se)$counts >= 1)
thresh_count = ceiling(0.1*(dim(filtered_se)[2]))
expressed_genes = names(which(expressed_counts > thresh_count))

#Filter genes based on biotype and chromosome
autosomes = c("1","10","11","12","13","14","15","16","17","18","19",
                      "2","20","21","22","3","4","5","6","7","8","9")
valid_gene_biotypes = c("lincRNA","protein_coding","IG_C_gene","IG_D_gene","IG_J_gene",
                        "IG_V_gene", "TR_C_gene","TR_D_gene","TR_J_gene", "TR_V_gene",
                        "3prime_overlapping_ncrna","known_ncrna", "processed_transcript",
                        "antisense","sense_intronic","sense_overlapping")
gene_data = rowData(filtered_se) %>% as.data.frame() %>% as_tibble()
valid_genes = dplyr::filter(gene_data, gene_type %in% valid_gene_biotypes, chromosome %in% autosomes)$gene_id

#Keep only valid and expressed genes
filtered_se = filtered_se[intersect(valid_genes, expressed_genes),]

#Normalise and make QTLtools matrix
cqn_se = normaliseSE_cqn(filtered_se)
qtltools_matrix = convertSEtoQTLtools(cqn_se, assay_name = "cqn")

#Export matrix to disk
write.table(qtltools_matrix, "processed/GEUVADIS/qtltools/phenotypes/GEUVADIS.featureCounts.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

#Extract sample list for the VCF file
columns = colnames(qtltools_matrix)
genotype_ids = columns[7:length(columns)]
write.table(genotype_ids, "processed/GEUVADIS/qtltools/phenotypes/GEUVADIS.featureCounts.samples", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


##### txrevise #####
#Export txrevise quants
txrevise_usage = readRDS("results/SummarizedExperiments/GEUVADIS_txrevise.rds")
txrevise_usage = txrevise_usage[rowData(txrevise_usage)$gene_id %in% rowData(cqn_se)$gene_id,]
filtered_usage = txrevise_usage[,txrevise_usage$qtl_mapping & txrevise_usage$superpopulation_code == "EUR"]

#Extract conditions and export to QTLtools
conditions = unique(filtered_usage$condition)
condition_list = setNames(as.list(conditions), conditions)
condition_se_list = purrr::map(condition_list, ~subsetSEByColumnValue(filtered_usage, "condition", .))

#Quantile normalise the ratios
qnorm_list = purrr::map(condition_se_list, ~normaliseSE_quantile(., assay_name = "usage"))

#Convert SE onbjects to QTLtools
qtltools_list = purrr::map(qnorm_list, ~convertSEtoQTLtools(., assay_name = "qnorm"))
saveQTLToolsMatrices(qtltools_list, output_dir = "processed/GEUVADIS/qtltools/input/txrevise/", file_suffix = "bed")

#Extract sample names
sample_names = purrr::map(qtltools_list, ~colnames(.)[-(1:6)])
saveQTLToolsMatrices(sample_names, output_dir = "processed/GEUVADIS/qtltools/input/txrevise/", file_suffix = "sample_names.txt", col_names = FALSE)

