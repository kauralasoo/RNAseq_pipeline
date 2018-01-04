library("dplyr")
rename_file = read.table("metadata/TwinsUK/rename_LCLs.txt", col.names = c("mv","old_bam","new_bam"), stringsAsFactors = FALSE)
sample_metadata = dplyr::select(rename_file, new_bam) %>% 
  tbl_df() %>% tidyr::separate(new_bam, c("sample_id", "suffix"), sep ="\\.") %>% 
  dplyr::select(sample_id) %>% 
  dplyr::mutate(fq1 = paste0(sample_id, "_1.fastq.gz"), fq2 = paste0(sample_id, "_2.fastq.gz"))

#Make snakemake strings
snakemake = dplyr::select(sample_metadata, sample_id, fq1, fq2) %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/TwinsUK/fastq/", fq1, ", processed/TwinsUK/fastq/", fq2,"]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakemake, "metadata/TwinsUK/TwinsUK_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
