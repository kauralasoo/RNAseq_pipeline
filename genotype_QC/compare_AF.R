geuvadis_var_info = importVariantInformation("coloc_data_temp/GEUVADIS_ge_LCL.variant_information.txt.gz")
kasela2017_var_info = importVariantInformation("coloc_data_temp/Kasela_2017_T-cell_CD4.variant_information.txt.gz")
kasela21 = dplyr::filter(kasela2017_var_info, chr == "21")

dplyr = dplyr::filter(kasela2017_var_info, chr == "21")
geuv21 = dplyr::filter(geuvadis_var_info, chr == "21")
multi_allelic = names(which(table(geuv21$pos) == 2))
geuv21 = dplyr::filter(geuv21, !(pos %in% multi_allelic))

shared_snps = dplyr::left_join(geuv21, kasela21, by = c("chr", "pos")) %>% dplyr::filter(!is.na(snp_id.y))
