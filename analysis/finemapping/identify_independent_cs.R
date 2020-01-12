library("purrr")
library("tidyr")
library("data.table")
library("devtools")
load_all("../eQTLUtils/")
library("igraph")

#Import metadata
HT12V4_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz",
                              col_types = "cccccdddcciidd")
ge_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz",
                          col_types = "cccccdddcciidd")
tx_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz",
                          col_types = "cccccdddcciidd")
txrev_meta = readr::read_tsv("https://zenodo.org/record/3366011/files/txrevise_Ensembl_96_phenotype_metadata.tsv.gz",
                             col_types = "cccccdddcciidd")
gene_ids = dplyr::bind_rows(ge_meta, tx_meta, txrev_meta, HT12V4_gene_ids) %>%
  dplyr::select(phenotype_id, gene_id, gene_name)

#Import credible sets
fairfax_2014_cs_df = importSusieCredibleSets("results/finemapping_results/Fairfax_2014/txt/") %>%
  dplyr::left_join(HT12V4_gene_ids, by = "phenotype_id")
blueprint = importSusieCredibleSets("results/finemapping_results/BLUEPRINT_SE/txt/") %>%
  dplyr::left_join(gene_ids, by = "phenotype_id")
alasoo = importSusieCredibleSets("results/finemapping_results/Alasoo_2018/txt/") %>%
  dplyr::left_join(gene_ids, by = "phenotype_id")

#Find independent credible sets
CD40_mono_cs = dplyr::filter(fairfax_2014_cs_df, qtl_group %like% "CL_0002057", gene_name == "CD40") %>%
  dplyr::mutate(cs_uid = paste(study, qtl_group, phenotype_id, cs_id, sep = "_")) %>%
  dplyr::group_by(cs_uid) %>%
  dplyr::mutate(cs_size = dplyr::n()) %>%
  dplyr::mutate(max_z = max(abs(z))) %>%
  dplyr::filter(max_z > 3.5)
ff_cs_list = setNames(dplyr::group_split(CD40_mono_cs), dplyr::group_keys(CD40_mono_cs)$cs_uid) %>%
  purrr::map(~.$variant_id)

blueprint_cs = dplyr::filter(blueprint, qtl_group %like% "monocyte", gene_id == "ENSG00000101017") %>%
  dplyr::mutate(cs_uid = paste(study, qtl_group, phenotype_id, cs_id, sep = "_")) %>%
  dplyr::group_by(cs_uid) %>%
  dplyr::mutate(cs_size = dplyr::n()) %>%
  dplyr::mutate(max_z = max(abs(z))) %>%
  dplyr::filter(max_z > 3.5)
bp_cs_list = setNames(dplyr::group_split(blueprint_cs), dplyr::group_keys(blueprint_cs)$cs_uid) %>%
  purrr::map(~.$variant_id)

alasoo_cs = dplyr::filter(alasoo, gene_id == "ENSG00000101017") %>%
  dplyr::mutate(cs_uid = paste(study, qtl_group, phenotype_id, cs_id, sep = "_")) %>%
  dplyr::group_by(cs_uid) %>%
  dplyr::mutate(cs_size = dplyr::n()) %>%
  dplyr::mutate(max_z = max(abs(z))) %>%
  dplyr::filter(max_z > 3.5)
alasoo_cs_list = setNames(dplyr::group_split(alasoo_cs), dplyr::group_keys(alasoo_cs)$cs_uid) %>%
  purrr::map(~.$variant_id)

#Split into list of credible sets
cs_list = c(ff_cs_list, bp_cs_list, alasoo_cs_list)
cs_metrics = calculatePairwiseCSMetrics(cs_list) %>%
  dplyr::filter(intersect > 0)

#Find connected compoents
self_pairs = dplyr::tibble(cs1 = names(cs_list), cs2 = names(cs_list))
overlap_pairs = dplyr::select(cs_metrics, cs1, cs2)
graph = igraph::graph_from_data_frame(dplyr::bind_rows(self_pairs, overlap_pairs), directed = F)
plot(graph)







#PTK2B example
PTK2B_mono_cs = dplyr::filter(cs_df_named, qtl_group %like% "CL_0002057", gene_name == "PTK2B") %>%
  dplyr::mutate(cs_uid = paste(study, qtl_group, phenotype_id, cs_id, sep = "_")) %>%
  dplyr::group_by(cs_uid) %>%
  dplyr::mutate(cs_size = dplyr::n()) %>%
  dplyr::mutate(max_z = max(abs(z))) %>%
  dplyr::filter(max_z > 3.5)

#Split into list of credible sets
cs_list = setNames(dplyr::group_split(PTK2B_mono_cs), dplyr::group_keys(PTK2B_mono_cs)$cs_uid) %>%
  purrr::map(~.$variant_id)
cs_metrics = calculatePairwiseCSMetrics(cs_list) %>%
  dplyr::filter(jaccard > 0.5)

#Find connected compoents
graph = igraph::graph_from_data_frame(cs_metrics, directed = F)
plot(graph)

