library(coloc)

eqtl_data = readRDS("coloc_data_temp/eqtl_data_raw_unfiltered.rds") %>%
  dplyr::filter(molecular_trait_id == "ENST00000466205") %>%
  dplyr::group_by(position) %>% 
  dplyr::mutate(alt_allele_count = length(alt)) %>% 
  dplyr::filter(alt_allele_count == 1)

gwas_data = readRDS("coloc_data_temp/gwas_data_raw_unfiltered.rds") %>%
  dplyr::filter(!is.na(ci_lower)) %>%
  dplyr::mutate(log_OR = log(odds_ratio)) %>%
  dplyr::mutate(se = (log(ci_upper)-log(ci_lower))/3.92)

shared_positions = intersect(eqtl_data$position, gwas_data$base_pair_location)

eqtl_shared = dplyr::filter(eqtl_data, position %in% shared_positions) %>% 
  dplyr::mutate(variant_id = as.character(position))
gwas_shared = dplyr::filter(gwas_data, base_pair_location %in% shared_positions) %>% 
  dplyr::mutate(variant_id = as.character(base_pair_location))


eQTL_dataset = list(pvalues = eqtl_shared$pvalue, 
                    N = 84, #The sample size of the eQTL dataset was 84
                    MAF = eqtl_shared$maf, 
                    type = "quant", 
                    beta = eqtl_shared$beta,
                    snp = eqtl_shared$variant_id)
gwas_dataset = list(beta = gwas_shared$log_OR, #If log_OR column is full of NAs then use beta column instead
                    varbeta = gwas_shared$se^2, 
                    type = "cc", 
                    snp = gwas_shared$variant_id,
                    s = 0.5, #This is acutally not used, because we already specified varbeta above.
                    MAF = eqtl_shared$maf)

coloc_res = coloc::coloc.abf(dataset1 = eQTL_dataset, dataset2 = gwas_dataset,p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
coloc_res$summary
