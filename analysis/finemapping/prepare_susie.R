library("dplyr")
library("devtools")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
load_all("../eQTLUtils/")


#Import qtltools input file:
qtlmap_out = "/gpfs/hpc/projects/eQTLCatalogue/qtlmap/RNAseq"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"


makeSusieInput <- function(qtltools_table, qtlmap_out, qcnorm_out){
  result = dplyr::mutate(qtltools_table, 
                         phenotype_list = file.path(qtlmap_out, "sumstats", paste0(qtl_subset, ".permuted.txt.gz")),
                         covariates = file.path(qtlmap_out, "PCA", qtl_subset, paste0(qtl_subset, ".covariates.txt"))) %>%
    dplyr::mutate(count_matrix = file.path(qcnorm_out, count_matrix)) %>%
    dplyr::transmute(qtl_subset, count_matrix, pheno_meta, 
                     sample_meta, vcf, phenotype_list, covariates)
}

#BLUEPRINT
qtltools_input = readr::read_tsv("results/finemap/BLUEPRINT_qtlmap_inputs.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/BLUEPRINT_study_file.tsv", sep = "\t", quote =F, row.names = F)


#Imputed studies
qtltools_input = readr::read_tsv("results/finemap/imputed_studies.txt", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/susie_study_file.tsv", sep = "\t", quote =F, row.names = F)

#Not imputed studies
qtltools_input = readr::read_tsv("results/finemap/not_imputed_studies.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/susie_study_file2.tsv", sep = "\t", quote =F, row.names = F)

#Array studies
qtltools_input = readr::read_tsv("results/finemap/array_input_files.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out = "/gpfs/hpc/projects/eQTLCatalogue/qtlmap/HumanHT12V4", qcnorm_out)
write.table(susie_input, "results/finemap/susie_study_file_microarray.tsv", sep = "\t", quote =F, row.names = F)

#TwinsUK
qtltools_input = readr::read_tsv("results/finemap/TwinsUK_qtlmap_inputs.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/susie_study_file_twins.tsv", sep = "\t", quote =F, row.names = F)

#CEDAR TopMed
qtltools_input = readr::read_tsv("results/finemap/CEDAR_2.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out = "/gpfs/hpc/projects/eQTLCatalogue/qtlmap/CEDAR_TOPmed/", qcnorm_out)
write.table(susie_input, "results/finemap/susie_CEDAR_TOPMed.tsv", sep = "\t", quote =F, row.names = F)


#Schimedel
qtlmap_out = "/gpfs/hpc/projects/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/results_130720"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"
qtltools_input = readr::read_tsv("results/finemap/Schmiedel_2018_qtlmap_inputs.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/Schmiedel_2018_finemap.tsv", sep = "\t", quote =F, row.names = F)

#Studies using GT field
qtlmap_out = "/gpfs/hpc/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"
qtltools_input = readr::read_tsv("results/finemap/GT_studies.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/GT_studies_finemap.tsv", sep = "\t", quote =F, row.names = F)

#Studies using DS field
qtlmap_out = "/gpfs/hpc/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"
qtltools_input = readr::read_tsv("results/finemap/DS_studies.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/DS_studies_finemap.tsv", sep = "\t", quote =F, row.names = F)

#GTEx ge only
qtlmap_out = "/gpfs/hpc/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/GTEx_v7"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"
qtltools_input = readr::read_tsv("results/finemap/ge_groups.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/ge_groups_finemap.tsv", sep = "\t", quote =F, row.names = F)


#Full GTEx V7 part 1
qtlmap_out = "/gpfs/hpc/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"
qtltools_input = readr::read_tsv("results/finemap/all_groups_part1.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/GTEx_part1.tsv", sep = "\t", quote =F, row.names = F)

#Full GTEx V7 part 2
qtlmap_out = "/gpfs/hpc/home/a72094/eQTLCatalogue/qtlmap/eQTL_Catalogue_r3/pipeline_out"
qcnorm_out = "/gpfs/hpc/projects/eQTLCatalogue/qcnorm/"
qtltools_input = readr::read_tsv("results/finemap/all_groups_part2.tsv", col_names = T)
susie_input = makeSusieInput(qtltools_input, qtlmap_out, qcnorm_out)
write.table(susie_input, "results/finemap/GTEx_part2.tsv", sep = "\t", quote =F, row.names = F)





