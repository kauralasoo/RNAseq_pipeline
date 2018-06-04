library("dplyr")

sample_meta = read.table("metadata/Alasoo_2018/RNA_sample_metadata.txt", stringsAsFactors = FALSE, header = TRUE)

snakakemake_str = dplyr::select(sample_meta, sample_id) %>% dplyr::tbl_df() %>% 
  dplyr::mutate(fq1 = paste0("processed/Macrophages_Alasoo_2018/fastq/", sample_id, ".1.fastq.gz"), 
                fq2 = paste0("processed/Macrophages_Alasoo_2018/fastq/", sample_id, ".2.fastq.gz")) %>%
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", fq1, ", ", fq2,"]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakakemake_str, "metadata/Alasoo_2018/Alasoo_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
