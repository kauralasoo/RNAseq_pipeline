
#Import gene metadata
array_gene_meta = read.table("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz", header = TRUE, stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
probe_gene_map = dplyr::select(array_gene_meta, phenotype_id, gene_name)
featureCounts_gene_meta = read.table("metadata/gene_metadata/featureCounts_Ensembl_92_gene_metadata.txt.gz", header = TRUE, stringsAsFactors = FALSE, sep = "\t") %>%
  dplyr::as_tibble()
gene_name_map = dplyr::select(featureCounts_gene_meta, phenotype_id, gene_name)


importSummaries <- function(summary_file, gene_names){
  extracted_results = read.table(summary_file, header = T, stringsAsFactors = F) %>%
    dplyr::as_tibble() %>%
    dplyr::select(qtl_group, phenotype_id, snp_id, p_nominal, beta, z_score, se) %>%
    dplyr::mutate(qtl_group = stringr::str_replace(qtl_group, "results/qtl_summary_stats/", "")) %>%
    tidyr::separate(qtl_group, into = c("study", "quant", "qtl_group"), sep = "\\/") %>%
    dplyr::left_join(gene_names, by = "phenotype_id")
  
  return(extracted_results)
}

#SLFN5 example
slfn5_summaries = importSummaries("results/extracted_variants/results.txt", gene_name_map)
slfn5 = dplyr::filter(extracted_results, phenotype_id == "ENSG00000166750")

ggplot(slfn5, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()


#featureCounts colocs
hits = importSummaries("results/extracted_variants/featureCounts_results.txt", gene_name_map)

#STARD10
gene = dplyr::filter(hits, gene_name == "STARD10", snp_id == "chr11_72759870_GGTTT_G")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

#TRAF1
gene = dplyr::filter(hits, gene_name == "TRAF1")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

#CARD9
#Colocalises in monocytes and LCLs, effect in the opposite direction. What is going on?
gene = dplyr::filter(hits, gene_name == "CARD9", snp_id == "chr9_136377295_G_A")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

#TNFRSF14 - exrtemely broad eQTL
gene = dplyr::filter(hits, gene_name == "TNFRSF14")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

#IKZF3
#Large effects in pancreatic islets and BLUEPRINT monocytes - caused by secondary QTLs?
gene = dplyr::filter(hits, gene_name == "IKZF3", snp_id == "chr17_39829548_A_G")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

#Microarray examples
hits_array = importSummaries("results/extracted_variants/array_results.txt", probe_gene_map)

#CTLA4
gene = dplyr::filter(hits_array, gene_name == "CTLA4")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()

#TRAF1
gene = dplyr::filter(hits_array, gene_name == "TRAF1")
ggplot(gene, aes(x = paste(qtl_group, study, sep = "-"), y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point(position=position_dodge(width=0.9)) + 
  geom_errorbar() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  coord_flip()



