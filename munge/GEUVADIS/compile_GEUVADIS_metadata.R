library("readr")

#Import metadata from ArrayExpress
metadata = read_tsv("metadata/GEUVADIS/GEAUVADIS_QC_passed_unique_mRNA.txt")
file_names = metadata[,c('Source Name', "Comment[FASTQ_URI]")]
colnames(file_names) = c("genotype_id", "fq")

#Import 1000G metadata
sample_meta = read_tsv("metadata/GEUVADIS/1000_genomes_sample_metadata.tsv", col_names =
           c("genotype_id", "sex", "biosample_id", "population_code", "population_name", 
             "superpopulation_code", "superpopulation_name", "data_collections"), skip = 1) %>%
  dplyr::select(-data_collections)

#Construct fq and sample names
geuvadis_metadata = dplyr::mutate(file_names, file_name = strsplit(fq, "/") %>% purrr::map_chr(., tail, 1) %>% unlist()) %>%
  tidyr::separate(file_name, c("sample_id", "suffix"), sep = "_", remove = FALSE) %>%
  dplyr::select(sample_id, genotype_id, file_name, suffix) %>% 
  dplyr::mutate(suffix = ifelse(suffix == "1.fastq.gz", "fq1", "fq2")) %>%
  tidyr::spread(suffix, file_name) %>%
  dplyr::left_join(sample_meta, by = "genotype_id")
write.table(geuvadis_metadata, "metadata/GEUVADIS/GEUVADIS_compiled_metadata.txt", sep="\t", quote=FALSE, row.names = FALSE)

#Make a snakemake string
snakemake = dplyr::select(geuvadis_metadata, sample_id, fq1, fq2) %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/GEUVADIS/fastq/", fq1, ", processed/GEUVADIS/fastq/", fq2,"]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakemake, "metadata/GEUVADIS/GEUVADIS_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
