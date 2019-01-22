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


