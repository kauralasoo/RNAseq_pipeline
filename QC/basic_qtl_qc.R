library("devtools")
load_all("../eQTLUtils/")

#### Alasoo_2018 ####

#Gene counts
naive_qtls = importQTLtoolsTable("results/Alasoo_2018/final/Alasoo_2018_ge_macrophage_naive/Alasoo_2018_ge_macrophage_naive.permuted.txt.gz")
dim(dplyr::filter(naive_qtls, p_fdr < 0.1))
hist(naive_qtls$p_perm)

#CD14
dplyr::filter(naive_qtls, group_id == "ENSG00000170458")$p_fdr

#SPOPL
dplyr::filter(naive_qtls, group_id == "ENSG00000144228")$p_fdr

#HMGCR
dplyr::filter(naive_qtls, group_id == "ENSG00000113161")$p_fdr


ifng_qtls = importQTLtoolsTable("results/Alasoo_2018/final/Alasoo_2018_ge_macrophage_IFNg/Alasoo_2018_ge_macrophage_IFNg.permuted.txt.gz")
dim(dplyr::filter(ifng_qtls, p_fdr < 0.1))
hist(ifng_qtls$p_perm)

#CD14
dplyr::filter(ifng_qtls, group_id == "ENSG00000170458")$p_fdr

#SPOPL
dplyr::filter(ifng_qtls, group_id == "ENSG00000144228")$p_fdr


#Exon counts
naive_exon_qtls = importQTLtoolsTable("results/Alasoo_2018/final/Alasoo_2018_exon_macrophage_naive/Alasoo_2018_exon_macrophage_naive.permuted.txt.gz")
dim(dplyr::filter(naive_exon_qtls, p_fdr < 0.1))
hist(naive_exon_qtls$p_perm)

#CD14
dplyr::filter(naive_exon_qtls, group_id == "ENSG00000170458")$p_fdr

#SPOPL
dplyr::filter(naive_exon_qtls, group_id == "ENSG00000144228")$p_fdr

#HMGCR
dplyr::filter(naive_exon_qtls, group_id == "ENSG00000113161")$p_fdr

#CD33
dplyr::filter(naive_exon_qtls, group_id == "ENSG00000105383")$p_fdr


### Transcript usage
naive_tu_qtls = importQTLtoolsTable("results/Alasoo_2018/final/Alasoo_2018_tx_macrophage_naive/Alasoo_2018_tx_macrophage_naive.permuted.txt.gz")
dim(dplyr::filter(naive_tu_qtls, p_fdr < 0.1))
hist(naive_tu_qtls$p_perm)

#CD14
dplyr::filter(naive_tu_qtls, group_id == "ENSG00000170458")$p_fdr

#HMGCR
dplyr::filter(naive_tu_qtls, group_id == "ENSG00000113161")$p_fdr

#CD33
dplyr::filter(naive_tu_qtls, group_id == "ENSG00000105383")$p_fdr

### Txrevise
naive_txrevise_qtls = importQTLtoolsTable("results/Alasoo_2018/final/Alasoo_2018_txrev_macrophage_naive/Alasoo_2018_txrev_macrophage_naive.permuted.txt.gz")
dim(dplyr::filter(naive_txrevise_qtls, p_fdr < 0.1))
hist(naive_txrevise_qtls$p_perm)

#HMGCR
dplyr::filter(naive_txrevise_qtls, group_id %like% "ENSG00000113161.contained")$p_fdr

dplyr::filter(naive_txrevise_qtls, group_id %like% "ENSG00000105383.contained")$p_fdr

#Identify genes with at least one QTL
naive_qtls = importQTLtoolsTable("results/finemapping/Alasoo_2018_permuted/final/Alasoo_2018_ge_macrophage_naive/Alasoo_2018_ge_macrophage_naive.permuted.txt.gz") %>%
  dplyr::filter(p_fdr < 0.1)
write.table(naive_qtls$phenotype_id, "results/finemapping/Alasoo_2018/naive_qtl_gene_list.txt", sep = "\t", quote = F, row.names = F, col.names = F)






