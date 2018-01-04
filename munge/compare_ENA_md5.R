library("dplyr")
library("tidyr")


#### GEUVADIS study ####
#Import ENA MD5 values
md5 = readr::read_tsv("metadata/GEUVADIS/ENA_md5.txt") %>%
  tidyr::separate(fastq_md5, c("md5_1","md5_2"), sep = ";") %>% 
  tidyr::separate(fastq_ftp, c("fq1","fq2"), sep = ";")

md5_df = dplyr::bind_rows(dplyr::transmute(md5, md5 = md5_1, fastq = fq1),
                          dplyr::transmute(md5, md5 = md5_2, fastq = fq2)) %>%
  dplyr::mutate(file_name = strsplit(fastq, "/") %>% purrr::map_chr(., tail, 1) %>% unlist())

#Import local MD5 values
md5_loc = read.table("metadata/GEUVADIS/GEUVADIS_md5.txt", stringsAsFactors = FALSE, col.names = c("md5", "file_name")) %>%
  tbl_df()
dplyr::anti_join(md5_df, md5_loc, by = c("file_name", "md5"))


#### CTCF study #####
md5 = readr::read_tsv("metadata/CTCF/CTCF_ENA_md5.txt") %>%
  tidyr::separate(fastq_md5, c("md5_1","md5_2"), sep = ";") %>% 
  tidyr::separate(fastq_ftp, c("fq1","fq2"), sep = ";")

md5_df = dplyr::bind_rows(dplyr::transmute(md5, md5 = md5_1, fastq = fq1),
                          dplyr::transmute(md5, md5 = md5_2, fastq = fq2)) %>%
  dplyr::mutate(file_name = strsplit(fastq, "/") %>% purrr::map_chr(., tail, 1) %>% unlist())

md5_loc = read.table("metadata/CTCF/CTCF_md5.txt", stringsAsFactors = FALSE, col.names = c("md5", "file_name")) %>%
  tbl_df()
dplyr::anti_join(md5_df, md5_loc, by = c("file_name", "md5"))


#### PU.1 study ####
md5 = readr::read_tsv("../Blood_ATAC/Blood_ATAC/data/Waszak_2015/E-MTAB-3657.sdrf.txt")
metadata = md5[,c(1,6,22,23,24)]
colnames(metadata) = c("sample_id", "genotype_id", "URL", "md5", "phenotype")

clean_metadata = dplyr::mutate(metadata, file_name = strsplit(URL, "/") %>% 
                                 purrr::map_chr(., tail, 1) %>% unlist()) %>%
  dplyr::select(sample_id, genotype_id, phenotype, file_name, md5, URL) %>%
  dplyr::mutate(phenotype = ifelse(phenotype == "PU.1", "PU1", phenotype)) %>% 
  dplyr::mutate(sample_id = paste(genotype_id, phenotype, sep = "_"))
write.table(clean_metadata, "../Blood_ATAC/Blood_ATAC/data/Waszak_2015/Waszak_2015_clean_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Extract PU.1 samples
pu1_data = dplyr::filter(clean_metadata, phenotype == "PU1")

#Compare md5 values
md5_loc = read.table("../Blood_ATAC/Blood_ATAC/data/Waszak_2015/PU1_md5.txt", stringsAsFactors = FALSE, col.names = c("md5", "file_name")) %>%
  tbl_df()
dplyr::anti_join(pu1_data, md5_loc, by = c("file_name", "md5"))

pu1_data

