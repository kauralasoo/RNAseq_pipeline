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
