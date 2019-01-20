library("oligo")
library("dplyr")
library("devtools")
library("stringr")
library("hugene10sttranscriptcluster.db")
load_all("../eQTLUtils/")

#Import gene metadata
transcript_meta = importBiomartMetadata("annotations/Ensembl92_biomart_download.txt.gz")
gene_meta = dplyr::select(transcript_meta, gene_id, chromosome, gene_start, gene_end, strand, gene_name, gene_type, gene_gc_content, gene_version) %>%     
  dplyr::mutate(phenotype_pos = as.integer(ceiling((gene_end + gene_start)/2))) %>%
  dplyr::distinct()
probe_metadata = as.data.frame(hugene10sttranscriptclusterENSEMBL) %>% as_tibble() %>%
  dplyr::transmute(phenotype_id = paste0("AFFY_", probe_id), quant_id = ensembl_id, group_id = ensembl_id, gene_id = ensembl_id) %>%
  dplyr::left_join(gene_meta, by = "gene_id") %>%
  dplyr::filter(chromosome %in% c(as.character(c(1:22)), "X","Y"))

#Keep only unique probe to gene mappings
unique_ids = dplyr::select(probe_metadata, phenotype_id) %>% 
  dplyr::group_by(phenotype_id) %>%
  dplyr::summarise(count = length(phenotype_id)) %>% 
  dplyr::arrange(-count) %>% 
  dplyr::filter(count == 1)
filtered_probes = dplyr::filter(probe_metadata, phenotype_id %in% unique_ids$phenotype_id)

#Export probe metadata to disk
gz1 = gzfile("metadata/gene_metadata/Affy_Human_Gene_1_0_ST_gene_metadata.txt.gz","w") 
write.table(filtered_probes, gz1, sep = "\t", quote = F, row.names = F)
close(gz1)

#Import and normalise CEL files
#cel_files = list.celfiles("results/ImmVar/cel_files/", listGzipped = TRUE) %>% file.path("results/ImmVar/cel_files/", .)
#affyRaw = read.celfiles(cel_files)
#eset <- rma(affyRaw)
#exp_matrix = exprs(eset)
#rownames(exp_matrix) = paste0("AFFY_", rownames(exp_matrix))
#exp_matrix = exp_matrix[probe_metadata$phenotype_id,]

#Import sample metadata
cd4 = readr::read_tsv("metadata/Raj_2014/ImmVar_CD4_sample_metadata.txt")
cd14 = readr::read_tsv("metadata/Raj_2014/ImmVar_CD14_sample_metadata.txt")
all_meta = dplyr::bind_rows(cd4, cd14)

mandatory_cols = c("sample_id", "genotype_id", "sex", "cell_type", "condition", "qtl_group", "timepoint", "read_length", "stranded", "paired", "protocol", "rna_qc_passed", "genotype_qc_passed","study")

useful_columns = c("!Sample_title","!Sample_geo_accession","!Sample_characteristics_ch1","!Sample_characteristics_ch1_1",
                   "!Sample_characteristics_ch1_2", "!Sample_characteristics_ch1_3", "!Sample_characteristics_ch1_4", 
                   "!Sample_characteristics_ch1_5")
selected_columns = all_meta[,useful_columns]
colnames(selected_columns) = c("study_sample_id", "sample_id", "age", "sex", "cell_type", "batch", "inclusion_markers", "exclusion_markers")

#Clean up the metdata
clean_meta = dplyr::mutate(selected_columns, age = stringr::str_replace(age,coll("age (yrs): "),"")) %>%
  dplyr::mutate(age = stringr::str_replace(age,coll("age: "),"")) %>%
  dplyr::mutate(age = as.integer(age)) %>%
  dplyr::mutate(sex = case_when(
    sex == "Sex: Female" ~ "female",
    sex == "gender: female" ~ "female",
    sex == "Sex: Male" ~ "male",
    sex == "gender: male" ~ "male",
  )) %>%
  dplyr::mutate(cell_type = case_when(
    cell_type == "cell type: monocytes from peripheral blood mononuclear cells (PBMCs)" ~ "monocyte",
    cell_type == "cell type: T4 Naive cells from human peripheral blood mononuclear cells (PBMCs)" ~ "T-cell"
   )) %>%
  dplyr::mutate(batch = stringr::str_replace(batch, fixed("batch: "), "batch_")) %>%
  dplyr::mutate(inclusion_markers = stringr::str_replace(inclusion_markers, fixed("inclusion markers: "), "")) %>%
  dplyr::mutate(exclusion_markers = stringr::str_replace(exclusion_markers, fixed("exclusion markers: "), "")) %>%
  dplyr::mutate(marker = case_when(
    cell_type == "T-cell" ~ "CD4",
    cell_type == "monocyte" ~ "CD14"
    )) %>%
  dplyr::mutate(qtl_group = paste(cell_type, marker, sep = "_")) %>%
  tidyr::separate(study_sample_id, c("genotype_id", "replicate"), "\\.", remove = FALSE) %>%
  dplyr::mutate(condition = "naive", timepoint = 0, read_length = NA, stranded = NA, paired = NA, 
                protocol = "hugene_10_ST", rna_qc_passed = TRUE, genotype_qc_passed = TRUE, 
                study = "Raj_2014") %>%
  dplyr::select(mandatory_cols, everything())
write.table(clean_meta, "metadata/cleaned/Raj_2014.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

#Import expression datq
eset = readRDS("results/ImmVar/ImmVar_hugene_eset.rds")
matrix = exprs(eset)
row.names(matrix) = paste0("AFFY_", row.names(matrix))
col_names = str_split(colnames(matrix), pattern = "\\_", n = 2, simplify = T)
colnames(matrix) = col_names[,1]

#Extract annotated subset
submatrix = matrix[filtered_probes$phenotype_id,clean_meta$sample_id]

#Export probe metadata to disk
gz1 = gzfile("results/expression_matrices/Human_Gene_1_0_ST/Raj_2014.tsv.gz","w") 
write.table(submatrix, gz1, sep = "\t", quote = F, row.names = T)
close(gz1)

#Make a Summarized experiment
se = makeSummarizedExperiment(assay = submatrix, row_data = filtered_probes, col_data = clean_meta, assay_name = "exprs")
raj_pca = transformSE_PCA(se, assay_name = "exprs", n_pcs = 10)


