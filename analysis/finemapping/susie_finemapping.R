library("GDSArray")
library("dplyr")
library("susieR")
library("purrr")

importVariantInformationFromGDS <- function(gdsfile){
  
  #Import individual columns
  snp_pos = GDSArray::GDSArray(gdsfile, "snp.position")
  snp_chromosome = GDSArray::GDSArray(gdsfile, "snp.chromosome")
  snp_id = GDSArray::GDSArray(gdsfile, "snp.rs.id")
  
  #Make a data frame
  snp_df = dplyr::data_frame(gds_snp_id = as.integer(names(snp_id)),
                             chromosome = as.vector(snp_chromosome),
                             pos = as.vector(snp_pos),
                             snp_id = as.vector(snp_id))
  return(snp_df)
}

extractGenotypeMatrixFromGDS <- function(chr, start, end, variant_information, gdsfile){
  
  #Extract variant ids from variant infromation
  var_meta = dplyr::filter(variant_information, chromosome == chr, pos > start, pos < end)
  gds_ids = var_meta$gds_snp_id
  
  #Extract genotype from the gds file
  geno = GDSArray::GDSArray(gdsfile, "genotype")
  genotype = as.matrix(geno[gds_ids,])
  rownames(genotype) = var_meta$snp_id
  
  return(genotype)
}

HumanHT_12_V4_gene_metadata = read.table("~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/HumanHT-12_V4_gene_metadata.txt.gz", sep = "\t", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
CEDAR = read.table("~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/CEDAR.tsv", stringsAsFactors = FALSE, header = TRUE, fill = TRUE)
CEDAR_normalized = read.table("~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/CEDAR.normalized.tsv.gz", sep = "\t", header = TRUE)
gdsfile = "~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/CEDAR_GRCh38.filtered.renamed.gds"
variant_info = importVariantInformationFromGDS(gdsfile)

#Import lead variants
monocytes = read.table("~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/eQTL_leads/monocyte_CD14.permuted.txt.gz", stringsAsFactors = FALSE, header = FALSE)
bcells = read.table("~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/eQTL_leads/B-cell_CD19.permuted.txt.gz", stringsAsFactors = FALSE, header = FALSE)
tcells8 = read.table("~/hpc/teaching/MTAT.03.239_Bioinformatics/projects/finemapping/eQTL_leads/T-cell_CD8.permuted.txt.gz", stringsAsFactors = FALSE, header = FALSE)

#Identify lead variants in each condition
mono_fdr = p.adjust(monocytes$V21, method = "fdr")
monocytes$FDR = mono_fdr
mono_lead_vars = dplyr::filter(monocytes, FDR < 0.01)

b_fdr = p.adjust(bcells$V21, method = "fdr")
bcells$FDR = b_fdr
b_lead_vars = dplyr::filter(bcells, FDR < 0.01)

t_fdr = p.adjust(tcells8$V21, method = "fdr")
tcells8$FDR = t_fdr
t_lead_vars = dplyr::filter(tcells8, FDR < 0.01)

newdf <- rbind(mono_lead_vars, b_lead_vars, t_lead_vars)
newdf = dplyr::arrange(newdf, FDR)
df = newdf[!duplicated(newdf$V6), ]
finemapping_probes = df$V6
fm_probes_list = setNames(as.list(finemapping_probes), finemapping_probes)

runFineMapping <- function(probe_id, HumanHT_12_V4_gene_metadata, CEDAR, CEDAR_normalized,
                           gdsfile, variant_info, ctype){
  
  print(probe_id)
  
  #Filter dataset by qtl_group
  filteredCedarByCellType = dplyr::filter(CEDAR,
                                          qtl_group == ctype,
                                          genotype_qc_passed == TRUE,
                                          rna_qc_passed == TRUE)
  sampleIds = filteredCedarByCellType$sample_id
  expressionsFilteredBySampleIds = CEDAR_normalized %>% select(one_of(sampleIds))
  
  #Perform finemapping for a single probe id
  geneData = dplyr::filter(HumanHT_12_V4_gene_metadata, phenotype_id == probe_id)
  position = geneData$phenotype_pos
    
  expressionInfos = expressionsFilteredBySampleIds[geneData$phenotype_id,]
  phenotype_data = t(expressionInfos)
  phenotype_mat = as.matrix(phenotype_data)
  colnames(phenotype_mat) = geneData$gene_name
  rownames(phenotype_mat) = filteredCedarByCellType$genotype_id
  
  phenotype_values = as.vector(phenotype_mat[,1])
  phenotype_values = (phenotype_values-mean(phenotype_values))/sd(phenotype_values)
  
  genotype_matrix = extractGenotypeMatrixFromGDS(chr = toString(geneData$chromosome),
                                                 start = position - 500000,
                                                 end = position + 500000,
                                                 variant_information = variant_info,
                                                 gdsfile = gdsfile)
  geno = genotype_matrix[,rownames(phenotype_mat)]
  geno_std = geno - apply(geno, 1, mean)
  geno_std = t(geno_std)
  fitted <- susie(geno_std, phenotype_values,
                  L = 10,
                  estimate_residual_variance = TRUE, 
                  estimate_prior_variance = FALSE,
                  scaled_prior_variance = 0.1,
                  verbose = TRUE)
  return(fitted)
}

# Wrap the finemapping function with purrr::safely so that an error in a 
# single finemapping call will not break the whole script
safeFineMapping = purrr::safely(runFineMapping)

#Finemap in three cell types
monocyte_results = purrr::map(fm_probes_list, ~safeFineMapping(.,HumanHT_12_V4_gene_metadata, 
                                                      CEDAR, CEDAR_normalized, gdsfile, variant_info, "monocyte_CD14"))
saveRDS(monocyte_results, "monocytes_cs.rds")
bcell_results = purrr::map(fm_probes_list, ~safeFineMapping(.,HumanHT_12_V4_gene_metadata, 
                                                     CEDAR, CEDAR_normalized, gdsfile, variant_info, "B-cell_CD19"))
saveRDS(bcell_results, "bcell_cs.rds")
tcell_results = purrr::map(fm_probes_list, ~safeFineMapping(.,HumanHT_12_V4_gene_metadata, 
                                             CEDAR, CEDAR_normalized, gdsfile, variant_info, "T-cell_CD8"))
saveRDS(tcell_results, "tcell_cs.rds")



