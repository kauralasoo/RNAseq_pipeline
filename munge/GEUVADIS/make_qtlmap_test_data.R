library("dplyr")
library("devtools")
library("SummarizedExperiment")
load_all("../eQTLUtils/")

#Import sample metadata
geuvadis = read.table("../SampleArcheology/studies/cleaned/GEUVADIS.tsv", 
               sep = "\t", stringsAsFactors = F, header =T) %>% dplyr::as_tibble()
ceu_meta = dplyr::filter(geuvadis, population_code == "CEU", rna_qc_passed, genotype_qc_passed)

#Import all lead variants
GEUVADIS_eQTLs = eQTLUtils::importQTLtoolsTable("results/qtl_summary_stats//GEUVADIS/featureCounts/LCL.permuted.txt.gz")

#Sample some eGenes and background genes
set.seed(1)
qtl_genes = dplyr::filter(GEUVADIS_eQTLs, p_nominal < 1e-12)$phenotype_id %>% sample(20)
eGenes = c(c("ENSG00000166750"), qtl_genes)

set.seed(1)
bg_genes = dplyr::filter(GEUVADIS_eQTLs, p_fdr > 0.5)$phenotype_id %>% sample(59)
selected_genes = c(eGenes, bg_genes)

#Import count matrix
GEUVADIS_counts = read.table("results/expression_matrices/featureCounts/GEUVADIS.tsv.gz")
geuvadis_gene_meta = read.table("metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz", stringsAsFactors = F, header = T) %>%
  dplyr::as_tibble()
#Normalize with cqn
geuvadis_se = makeSummarizedExperiment(GEUVADIS_counts, geuvadis_gene_meta, geuvadis, "counts")
geuvadis_cqn = normaliseSE_cqn(geuvadis_se)

#Subset by samples and genes
cqn_subset = geuvadis_cqn[selected_genes, ceu_meta$sample_id]
subset_mat = assays(cqn_subset)$cqn
subset_df = dplyr::mutate(as.data.frame(subset_mat), phenotype_id = rownames(subset_mat)) %>%
  dplyr::as_tibble() %>% dplyr::select(phenotype_id, everything())
write.table(subset_df, "results/qtlmap_test/GEUVADIS_cqn.tsv", sep = "\t", quote = F, row.names = F)

#Extract minimal sample metadata
min_sample_meta = colData(cqn_subset) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  dplyr::select(sample_id, genotype_id, qtl_group)
write.table(min_sample_meta, "results/qtlmap_test/GEUVADIS_sample_metadata.tsv", sep = "\t", quote = F, row.names = F)

#Extract phenotype metadata
min_phenotype_meta = rowData(cqn_subset) %>% 
  as.data.frame() %>% 
  as_tibble() %>% 
  dplyr::select(phenotype_id, group_id, chromosome, phenotype_pos, strand)
write.table(min_phenotype_meta, "results/qtlmap_test/GEUVADIS_phenotype_metadata.tsv", sep = "\t", quote = F, row.names = F)

#Make regions file to extract genotypes from VCF
regions = dplyr::select(min_phenotype_meta, chromosome, phenotype_pos) %>%
  dplyr::mutate(start = phenotype_pos - 100000, end = phenotype_pos + 100000) %>%
  dplyr::select(chromosome, start, end) %>%
  dplyr::arrange(chromosome, start) %>%
  dplyr::mutate(start = ifelse(start < 0, 0, start))
write.table(regions, "results/qtlmap_test/GEUVADIS_genotype_regions.txt", sep = "\t", quote = F, row.names = F, col.names = F)
  



