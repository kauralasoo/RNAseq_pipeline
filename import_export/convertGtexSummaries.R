library("data.table")

#Import data
tpms = readr::read_tsv("results/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct.gz", skip = 2)
attr = readr::read_tsv("results/GTEx/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt")
tissue_meta = dplyr::select(attr, SMTSD, SMUBRID) %>% 
  dplyr::distinct() %>%
  dplyr::transmute(gtex_tissue = SMTSD, tissue_ontology_id = ifelse(SMUBRID %like% "EFO", SMUBRID, paste0("UBER_", SMUBRID)))
  
