library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cqn")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
library("devtools")
library("stringr")
library("data.table")
load_all("../eQTLUtils/")

#Specify mandatory metadata columns
mandatory_cols = c("sample_id", "genotype_id", "sex", "cell_type", "condition", "qtl_group", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")

#### Fairfax_2014 ####
#Import SE
fairfax_2014 = readRDS("results/SummarizedExperiments/microarray/Fairfax_2014.rds")

#Import Plink sex estimates
set1 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set1.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
set2 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set2.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
plink_sex = dplyr::bind_rows(set1, set2) %>% 
  dplyr::select(genotype_id, plink_sex) %>% 
  dplyr::mutate(sex = "NA") %>% 
  dplyr::mutate(sex = ifelse(plink_sex == 1, "male", sex)) %>%
  dplyr::mutate(sex = ifelse(plink_sex == 2, "female", sex)) %>%
  dplyr::select(genotype_id, sex)

#Reformat sample metadata
sample_metadata = colData(fairfax_2014) %>% as.data.frame() %>% as_tibble() %>%
  dplyr::select(-sex) %>%
  dplyr::left_join(plink_sex, by = "genotype_id") %>%
  dplyr::rename(genotype_qc_passed = genotype_QC_passed) %>%
  dplyr::rename(rna_qc_passed = RNA_QC_passed) %>%
  dplyr::mutate(genotype_qc_passed = ifelse(is.na(sex), FALSE, TRUE)) %>%
  dplyr::mutate(genotype_id = paste0("FF14_", genotype_id)) %>%
  dplyr::mutate(condition = ifelse(condition == "naive", condition, paste0(condition, timepoint))) %>%
  dplyr::mutate(cell_type = "monocyte") %>%
  dplyr::mutate(qtl_group = paste(cell_type, condition, sep = "_")) %>%
  dplyr::mutate(protocol = "HumanHT-12_V4") %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/Fairfax_2014.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Extract gene metadata
gene_metadata = rowData(fairfax_2014) %>% as.data.frame() %>%
  dplyr::mutate(quant_id = gene_id, group_id = gene_id)
gz1 = gzfile("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz","w") 
write.table(gene_metadata, gz1, sep = "\t", quote = F, row.names = F)
close(gz1)

#Save expression matrix
mat = assays(fairfax_2014)$exprs
gz2 = gzfile("results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz1)


#### Fairfax_2012 ####
#Import SE
fairfax_2012 = readRDS("results/SummarizedExperiments/microarray/Fairfax_2012.rds")

#Import Plink sex estimates
set1 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set1.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
set2 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set2.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
plink_sex = dplyr::bind_rows(set1, set2) %>% 
  dplyr::select(genotype_id, plink_sex) %>% 
  dplyr::mutate(sex = "NA") %>% 
  dplyr::mutate(sex = ifelse(plink_sex == 1, "male", sex)) %>%
  dplyr::mutate(sex = ifelse(plink_sex == 2, "female", sex)) %>%
  dplyr::select(genotype_id, sex)

#Reformat 
sample_metadata = colData(fairfax_2012) %>% as.data.frame() %>% as_tibble() %>%
  dplyr::select(-sex) %>% 
  dplyr::mutate(genotype_id = as.integer(stringr::str_replace(genotype_id, "b_", ""))) %>%
  dplyr::left_join(plink_sex, by = "genotype_id") %>%
  dplyr::mutate(qtl_group = paste(cell_type, marker, sep = "_")) %>%
  dplyr::rename(genotype_qc_passed = genotype_QC_passed) %>%
  dplyr::rename(rna_qc_passed = RNA_QC_passed) %>%
  dplyr::mutate(genotype_id = paste0("FF14_", genotype_id)) %>%
  dplyr::mutate(protocol = "HumanHT-12_V4") %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/Fairfax_2012.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Save expression matrix
mat = assays(fairfax_2012)$exprs
gz2 = gzfile("results/expression_matrices/HumanHT-12_V4/Fairfax_2012.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)



#### Naranbhai_2012 ####
naranbhai_2015 = readRDS("results/SummarizedExperiments/microarray/Naranbhai_2015.rds")

#Import Plink sex estimates
set1 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set1.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
set2 = readr::read_delim("metadata/Fairfax_2014/Fairfax_2014_set2.fam", delim = " ", col_names = c("family_id", "genotype_id", "father", "mother", "plink_sex", "phenotype"))
plink_sex = dplyr::bind_rows(set1, set2) %>% 
  dplyr::select(genotype_id, plink_sex) %>% 
  dplyr::mutate(sex = "NA") %>% 
  dplyr::mutate(sex = ifelse(plink_sex == 1, "male", sex)) %>%
  dplyr::mutate(sex = ifelse(plink_sex == 2, "female", sex)) %>%
  dplyr::select(genotype_id, sex)

#Reformat metadata
sample_metadata = colData(naranbhai_2015) %>% as.data.frame() %>% as_tibble() %>%
  dplyr::mutate(sample_id = stringr::str_replace(sample_id, "Sample ", "CD16_")) %>%
  dplyr::select(-sex) %>% 
  dplyr::mutate(cell_type = "neutrophil") %>%
  dplyr::mutate(genotype_id = as.integer(genotype_id)) %>%
  dplyr::left_join(plink_sex, by = "genotype_id") %>%
  dplyr::mutate(genotype_id = paste0("FF14_", genotype_id)) %>%
  dplyr::mutate(qtl_group = paste(cell_type, marker, sep = "_")) %>%
  dplyr::rename(genotype_qc_passed = genotype_QC_passed) %>%
  dplyr::rename(rna_qc_passed = RNA_QC_passed) %>%
  dplyr::mutate(protocol = "HumanHT-12_V4") %>%
  dplyr::select(mandatory_cols, everything())
  
write.table(sample_metadata, "metadata/cleaned/Naranbhai_2015.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Save expression matrix
mat = assays(naranbhai_2015)$exprs
colnames(mat) = sample_metadata$sample_id
gz2 = gzfile("results/expression_matrices/HumanHT-12_V4/Naranbhai_2012.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)


#### CEDAR ####
cedar_se = readRDS("results/SummarizedExperiments/microarray/CEDAR.rds")

#Perform cedar QC for non-blood tissues
cedar_filtered = cedar_se[, cedar_se$cell_type %in% c("rectum", "ileum", "transverse colon", "monocytes")]

#Normalize and regress out batch effects
cedar_norm = array_normaliseSE(cedar_filtered, norm_method = "quantile", assay_name = "exprs", log_transform = TRUE, adjust_batch = TRUE, filter_quality = FALSE)

#Perform PCA
pca = eQTLUtils::transformSE_PCA(cedar_norm, assay_name = "norm_exprs")
ggplot(pca$pca_matrix, aes(x = PC1, y = PC2, color = cell_type, label = sample_id)) + geom_text()
ggplot(dplyr::filter(pca$pca_matrix, cell_type == "ileum"), aes(x = PC1, y = PC2, color = cell_type, label = sample_id)) + geom_text()
ggplot(dplyr::filter(pca$pca_matrix, cell_type == "rectum"), aes(x = PC1, y = PC2, color = cell_type, label = sample_id)) +  geom_point()
ggplot(dplyr::filter(pca$pca_matrix, cell_type == "transverse colon"), aes(x = PC1, y = PC2, color = cell_type, label = sample_id)) + geom_text()

pca_matrix = dplyr::filter(pca$pca_matrix, !(sample_id %in% outlier_samples))
ggplot(pca_matrix, aes(x = PC3, y = PC4, color = cell_type, label = sample_id)) + geom_point()


#Import genotype metadata
genotype_meta = readr::read_delim("metadata/CEDAR/CEDAR_genotype_metadata.sdrf.txt", delim = "\t")
genotype_meta = genotype_meta[,c("Source Name","Characteristics[individual]","Characteristics[age]","Characteristics[clinical history]","Characteristics[excluded_ind]")]
colnames(genotype_meta) = c("genotype_array_id", "genotype_id", "age", "smoker", "included") 

#Import names from the VCF file
genotyped = read.table("metadata/CEDAR/CEDAR_genotyped_names.txt")

#Filter metadata
genotype_meta = dplyr::mutate(genotype_meta, included = ifelse(included == "excluded from study", FALSE, TRUE)) %>%
  dplyr::distinct() %>%
  dplyr::filter(genotype_array_id %in% genotyped$V1)

#Export genotype-name map
map = dplyr::select(genotype_meta, genotype_array_id, genotype_id)
write.table(map, "metadata/CEDAR/CEDAR_genotype_name_map.txt", sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Define outlier samples
rec_samples = c("REIPC108","REIPC067","REIPC060","RE1IPC260_bis","RE1IPC434_bis","RE3IPC278_bis","RE1IPC434","RE3IPC361",
                "RE3IPC278","RE3IPC263","RE3IPC260","RE3IPC289","RE3IPC262","REIPC257","RE3IPC253")
il_samples = c("IL1IPC246","IL1IPC375","IL1IPC371","IL1IPC306","IL1IPC305","IL1IPC434","IL1IPC427","IL1IPC312","IL1IPC309",
               "IL1IPC370","IL1IPC369","IL1IPC301","IL1IPC300","IL1IPC368","IL1IPC353","IL1IPC298","IL1IPC297","IL1IPC350",
               "IL1IPC348","IL1IPC293","IL1IPC289","IL1IPC347","IL1IPC344","IL1IPC288","IL1IPC282","IL1IPC368_bis",
               "IL2IPC106", "ILIPC152", "ILIPC152bis", "IL1IPC315", "IL1IPC326")
colon_samples = c("TRIPC169","TR1IPC427","TR1IPC280","TR1IPC399","TR1IPC278","TR1IPC368")
outlier_samples = c(rec_samples,il_samples,colon_samples)

#Extract sample metadata
sample_metadata = colData(cedar_se) %>% as.data.frame() %>% as_tibble() %>%
  dplyr::mutate(marker = ifelse(marker %like% "CD", marker, "none")) %>%
  dplyr::mutate(cell_type = case_when(
    cell_type == "monocytes" ~ "monocyte",
    cell_type == "neutrophils" ~ "neutrophil",
    cell_type == "transverse colon" ~ "transverse_colon",
    TRUE ~ cell_type
  )) %>%
  dplyr::mutate(qtl_group = ifelse(marker == "none", cell_type, paste(cell_type, marker, sep = "_"))) %>%
  dplyr::mutate(protocol = "HumanHT-12_V4") %>%
  dplyr::rename(genotype_qc_passed = genotype_QC_passed) %>%
  dplyr::rename(rna_qc_passed = RNA_QC_passed) %>%
  dplyr::left_join(genotype_meta, by = "genotype_id") %>%
  dplyr::mutate(rna_qc_passed = ifelse(cell_type %in% c("ileum","rectum", "transverse_colon"), TRUE, rna_qc_passed)) %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %in% outlier_samples, FALSE, rna_qc_passed)) %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %like% "bis", FALSE, rna_qc_passed)) %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %like% "quat", FALSE, rna_qc_passed)) %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %like% "sex", FALSE, rna_qc_passed)) %>%
  dplyr::mutate(rna_qc_passed = ifelse(sample_id %like% "quin", FALSE, rna_qc_passed)) %>%
  dplyr::mutate(included = ifelse(is.na(included), FALSE, included)) %>%
  dplyr::mutate(genotype_qc_passed = ifelse(included, TRUE, FALSE)) %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/CEDAR.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Save expression matrix
mat = round(assays(cedar_se)$exprs, 1)
gz2 = gzfile("results/expression_matrices/HumanHT-12_V4/CEDAR.tsv.gz", "w")
write.table(mat, gz2, sep = "\t", quote = FALSE)
close(gz2)






#### Kasela_2017 ####
kasela_se = readRDS("results/SummarizedExperiments/microarray/Kasela_2017.rds")

sample_metadata = colData(kasela_se) %>% as.data.frame() %>% as_tibble() %>%
  dplyr::rename(genotype_qc_passed = genotype_QC_passed) %>%
  dplyr::rename(rna_qc_passed = RNA_QC_passed) %>%
  dplyr::mutate(qtl_group = paste(cell_type, marker, sep = "_")) %>%
  dplyr::mutate(protocol = "HumanHT-12_V4") %>%
  dplyr::select(mandatory_cols, everything())
write.table(sample_metadata, "metadata/cleaned/Kasela_2017.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

  

