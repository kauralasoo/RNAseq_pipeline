library("data.table")
library("devtools")
load_all("../eQTLUtils/")

#Import gene metadata
gene_metadata = read.table("~/annotations/eQTLCatalogue/v0.1/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz", stringsAsFactors = F, header = T) %>%
  dplyr::as.tbl() 
#Specify required phenotype metadata columns
required_phenotype_meta_columns = c("phenotype_id","quant_id","group_id","gene_id","chromosome","gene_start",
                                    "gene_end","strand","gene_name","gene_type","gene_version","phenotype_pos")
required_gene_meta_columns = c(required_phenotype_meta_columns, "phenotype_gc_content", "phenotype_length")


### txrevise events ####
list = as.list(c("results/TPM/txrevise.grp_1.contained.TPM.merged.txt", 
                 "results/TPM/txrevise.grp_2.contained.TPM.merged.txt",
                 "results/TPM/txrevise.grp_1.upstream.TPM.merged.txt",
                 "results/TPM/txrevise.grp_2.upstream.TPM.merged.txt",
                 "results/TPM/txrevise.grp_1.downstream.TPM.merged.txt",
                 "results/TPM/txrevise.grp_2.downstream.TPM.merged.txt"))
event_quants = purrr::map_df(list, ~readr::read_tsv(.))

#Extract relevant infromation from the txrevise events
txrevise_meta = dplyr::select(event_quants, phenotype_id) %>% 
  tidyr::separate(phenotype_id, c("gene_id", "grp", "pos", "transcript"), sep = "\\.", remove = F) %>%
  dplyr::mutate(quant_id = paste(gene_id, grp, pos, sep = "."), group_id = paste(gene_id, pos, sep = ".")) %>% 
  dplyr::select(phenotype_id, quant_id, group_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -phenotype_gc_content, -group_id, -quant_id, -phenotype_length), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("results/txrevise_Ensembl_96_CAGE_10bp_phenotype_metadata.tsv.gz", "w")
write.table(txrevise_meta, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

