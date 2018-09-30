
EGA_idat = read.table("metadata/HipSci/HipSci_EGA_idat_files.txt")
idat = read.table("metadata/HipSci/idat_files.txt")

#Find files that failed to download
dplyr::anti_join(EGA_idat, idat, by = "V1")


idat_meta = tidyr::separate(idat, V1, c("genotype_id", "genotype_array","chip_id", "dummy", "date", "dummy2"), sep = "\\.") %>% 
  dplyr::as_tibble() %>% 
  dplyr::select(-dummy, -dummy2, genotype_array)

missing_files = dplyr::group_by(idat_meta, genotype_id) %>% dplyr::mutate(file_count = length(genotype_id)) %>% dplyr::filter(file_count < 2)
