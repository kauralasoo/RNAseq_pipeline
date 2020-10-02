library("dplyr")
library("gwasvcf")

eqtl_data = readr::read_tsv("results/figure_data/BLUEPRINT_PE.T-cell_ge.nominal.sorted.tsv.gz") %>%
  dplyr::filter(molecular_trait_id == "ENSG00000153250")
gwas_data = gwasvcf::query_gwas("results/figure_data/LC_GWAS_subset.GRCh38.sorted.vcf.gz", chrompos = "2:160100000-160700000") %>% 
  gwasvcf::vcf_to_granges() %>% 
  as.data.frame() %>%
  dplyr::as_tibble() %>%
  dplyr::mutate(position = start)
gwas_data_selected = dplyr::semi_join(gwas_data, eqtl_data, by = "position")

#Import eQTL credible sets
cs = readr::read_tsv("results/figure_data/BLUEPRINT_PE.T-cell_ge.purity_filtered.txt.gz") %>%
  dplyr::filter(phenotype_id == "ENSG00000153250")

#Flag the credible set in GWAS results
gwas_flagged = dplyr::mutate(gwas_data_selected, in_cs = ifelse(position %in% cs$pos, TRUE, FALSE))

#Make manhattan plots
ggplot(eqtl_data, aes(y = -log(pvalue), x = position)) + geom_point()
ggplot(gwas_flagged, aes(y = LP, x = position, color = in_cs)) + geom_point()

#Import all coloc results
files_list = list.files("results/figure_data/LC-ebi-a-GCST004627/")
paths = paste0("results/figure_data/LC-ebi-a-GCST004627/", files_list)
path_list = setNames(as.list(paths), files_list)
coloc_res = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "dataset")
rbms1_colocs = dplyr::filter(coloc_res, molecular_trait_id == "ENSG00000153250", variant == "chr2_160468964_A_T")


#Import all coloc results
files_list = list.files("results/figure_data/LC-ebi-a-GCST004627_GTEx_V8/")
paths = paste0("results/figure_data/LC-ebi-a-GCST004627_GTEx_V8/", files_list)
path_list = setNames(as.list(paths), files_list)
coloc_res = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "dataset")
rbms1_colocs = dplyr::filter(coloc_res, molecular_trait_id == "ENSG00000153250", variant == "chr2_160468964_A_T")




