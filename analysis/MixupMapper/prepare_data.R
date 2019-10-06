suppressPackageStartupMessages(library("devtools"))
load_all("../eQTLUtils/")
library("SummarizedExperiment")


extractExpressionMatrix <- function(se, qtl_qroup){
  selected_se = se[,se$qtl_group == qtl_qroup]
  matrix = assays(selected_se)$counts
  df = matrix %>% as.data.frame() %>%
    dplyr::as_tibble() %>%
    dplyr::mutate(phenotype_id = rownames(matrix)) %>%
    dplyr::select(phenotype_id, everything()) 
  return(df)
}

makeCouplings <- function(se, qtl_group){
  selected_se = se[,se$qtl_group == qtl_group]
  res = colData(selected_se) %>%
    as.data.frame() %>%
    as_tibble() %>%
    dplyr::select(genotype_id, sample_id)
  return(res)
}

opt = list(
  expression_matrix = "~/datasets/processed/HumanHT-12_V4/CEDAR.HumanHT-12_V4_norm_exprs.tsv.gz",
  sample_meta = "~/projects/SampleArcheology/studies/cleaned/CEDAR.tsv",
  phenotype_meta = "~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz"
)

#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                        row_data = phenotype_meta, 
                                                        col_data = sample_metadata, 
                                                        quant_method = "gene_counts",
                                                        reformat = FALSE)

#Extraxt all qtl_groups
qtl_groups = unique(se$qtl_group)
group_list = setNames(as.list(qtl_groups), qtl_groups)

#Extact matrices
matrices = purrr::map(group_list, ~extractExpressionMatrix(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".tsv"))), qtl_groups)
purrr::map2(matrices, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Make genotype-phenotype couplings
couplings = purrr::map(group_list, ~makeCouplings(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".geno_pheno_coupling.tsv"))), qtl_groups)
purrr::map2(couplings, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Make phenotype metada
res = rowData(monocyte_se) %>%
  as.data.frame() %>%
  as_tibble() %>%
  dplyr::transmute(Platform = "HT12v4", HT12v4_ArrayAddress = phenotype_id, Symbol = gene_name,
                   Chr = chromosome, ChrStart = phenotype_pos, ChrEnd = phenotype_pos+1, Probe = phenotype_id,
                   Seq = gene_id)
write.table(res, "results/mixupmapper/phenotype_metadata.tsv", sep = "\t", quote = F, row.names = F)

#Import all lead variants
files = list.files("results/lead_variants/CEDAR/", full.names = T)
files_list = setNames(as.list(files), files)
leads = purrr::map_df(files_list, ~eQTLUtils::importQTLtoolsNominalTable(.)) 
snp_list = dplyr::filter(leads, p_fdr < 0.05) %>% dplyr::select(snp_id) %>% dplyr::distinct()
write.table(snp_list, "results/mixupmapper/snp_list.txt", sep = "\t", quote =F, row.names = F, col.names = F)

#Import neutrophil leads from Naranbhai_2015
leads = eQTLUtils::importQTLtoolsNominalTable("results/lead_variants/Naranbhai_2015/neutrophil_CD16.lead_variants.txt") 
snp_list = dplyr::filter(leads, p_fdr < 0.05) %>% dplyr::select(snp_id) %>% dplyr::distinct()
write.table(snp_list, "results/mixupmapper/neutrophil_snp_list.txt", sep = "\t", quote =F, row.names = F, col.names = F)


#Naranbhai_2015
opt = list(
  expression_matrix = "~/datasets/processed/HumanHT-12_V4/Naranbhai_2015.HumanHT-12_V4_norm_exprs.tsv",
  sample_meta = "~/projects/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv",
  phenotype_meta = "~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz"
)

#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                        row_data = phenotype_meta, 
                                                        col_data = sample_metadata, 
                                                        quant_method = "gene_counts",
                                                        reformat = FALSE)

#Extraxt all qtl_groups
qtl_groups = unique(se$qtl_group)
group_list = setNames(as.list(qtl_groups), qtl_groups)

#Extact matrices
matrices = purrr::map(group_list, ~extractExpressionMatrix(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".tsv"))), qtl_groups)
purrr::map2(matrices, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Make genotype-phenotype couplings
couplings = purrr::map(group_list, ~makeCouplings(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".geno_pheno_coupling.tsv"))), qtl_groups)
purrr::map2(couplings, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Fairfax_2012
opt = list(
  expression_matrix = "~/datasets/processed/HumanHT-12_V4/Fairfax_2012.HumanHT-12_V4_norm_exprs.tsv",
  sample_meta = "~/projects/SampleArcheology/studies/cleaned/Fairfax_2012.tsv",
  phenotype_meta = "~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz"
)

#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                        row_data = phenotype_meta, 
                                                        col_data = sample_metadata, 
                                                        quant_method = "gene_counts",
                                                        reformat = FALSE)

#Extraxt all qtl_groups
qtl_groups = unique(se$qtl_group)
group_list = setNames(as.list(qtl_groups), qtl_groups)

#Extact matrices
matrices = purrr::map(group_list, ~extractExpressionMatrix(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".tsv"))), qtl_groups)
purrr::map2(matrices, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Make genotype-phenotype couplings
couplings = purrr::map(group_list, ~makeCouplings(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".geno_pheno_coupling.tsv"))), qtl_groups)
purrr::map2(couplings, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Fairfax_2014
opt = list(
  expression_matrix = "~/datasets/processed/HumanHT-12_V4/Fairfax_2014.HumanHT-12_V4_norm_exprs.tsv",
  sample_meta = "~/projects/SampleArcheology/studies/cleaned/Fairfax_2014.tsv",
  phenotype_meta = "~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz"
)

#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                        row_data = phenotype_meta, 
                                                        col_data = sample_metadata, 
                                                        quant_method = "gene_counts",
                                                        reformat = FALSE)

#Extraxt all qtl_groups
qtl_groups = unique(se$qtl_group)
group_list = setNames(as.list(qtl_groups), qtl_groups)

#Extact matrices
matrices = purrr::map(group_list, ~extractExpressionMatrix(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".tsv"))), qtl_groups)
purrr::map2(matrices, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Make genotype-phenotype couplings
couplings = purrr::map(group_list, ~makeCouplings(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".geno_pheno_coupling.tsv"))), qtl_groups)
purrr::map2(couplings, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Kasela_2017
opt = list(
  expression_matrix = "~/datasets/processed/HumanHT-12_V4/Kasela_2017.HumanHT-12_V4_norm_exprs.tsv",
  sample_meta = "~/projects/SampleArcheology/studies/cleaned/Kasela_2017.tsv",
  phenotype_meta = "~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz"
)

#Import all files
expression_matrix = readr::read_tsv(opt$expression_matrix)
sample_metadata = utils::read.csv(opt$sample_meta, sep = '\t', stringsAsFactors = F)
phenotype_meta = readr::read_delim(opt$phenotype_meta, delim = "\t", col_types = "ccccciiicciidi")

#Make a SummarizedExperiment of the expression data
se = eQTLUtils::makeSummarizedExperimentFromCountMatrix(assay = expression_matrix, 
                                                        row_data = phenotype_meta, 
                                                        col_data = sample_metadata, 
                                                        quant_method = "gene_counts",
                                                        reformat = FALSE)

#Extraxt all qtl_groups
qtl_groups = unique(se$qtl_group)
group_list = setNames(as.list(qtl_groups), qtl_groups)

#Extact matrices
matrices = purrr::map(group_list, ~extractExpressionMatrix(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".tsv"))), qtl_groups)
purrr::map2(matrices, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))

#Make genotype-phenotype couplings
couplings = purrr::map(group_list, ~makeCouplings(se,.))
path_list = setNames(as.list(file.path("results/mixupmapper", paste0(qtl_groups,".geno_pheno_coupling.tsv"))), qtl_groups)
purrr::map2(couplings, path_list, ~write.table(.x, .y, sep = "\t", quote = FALSE, row.names = FALSE))


