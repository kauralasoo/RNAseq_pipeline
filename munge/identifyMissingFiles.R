library("dplyr")

#GEUVADIS
#Import wget command
files = readr::read_delim("metadata/GEUVADIS/GEUVADIS_download_fastq.sh", col_names = c("wget","URL"), delim = " ")
downloaded_files = readr::read_delim("metadata/GEUVADIS/GEUVADIS_downloaded.txt", delim = " ", col_name = c("file_name"))

#Split into file names
file_names = dplyr::mutate(files, file_name = strsplit(URL, "/") %>% purrr::map_chr(., tail, 1) %>% unlist())

#Find missing files
missing_files = dplyr::anti_join(file_names, downloaded_files, by = "file_name") %>% dplyr::select(-file_name)
write.table(missing_files, "metadata/GEUVADIS/GEUVADIS_missing.sh", sep = " ", col.names=FALSE, row.names = FALSE, quote = FALSE)


#Nedelec
files = readr::read_delim("metadata/Nedelec_2017/Macrophages_Nedelec_2017.sh", col_names = c("wget","URL"), delim = " ")
downloaded_files = readr::read_delim("metadata/Nedelec_2017/Nedelec_downloaded.txt", delim = " ", col_name = c("file_name"))

#Split into file names
file_names = dplyr::mutate(files, file_name = strsplit(URL, "/") %>% purrr::map_chr(., tail, 1) %>% unlist())

#Find missing files
missing_files = dplyr::anti_join(file_names, downloaded_files, by = "file_name") %>% dplyr::select(-file_name)
write.table(missing_files, "metadata/GEUVADIS/GEUVADIS_missing.sh", sep = " ", col.names=FALSE, row.names = FALSE, quote = FALSE)


#PU.1
files = readr::read_delim("../Blood_ATAC/Blood_ATAC/download_scripts/PU1.sh", col_names = c("wget","URL"), delim = " ")
downloaded_files = readr::read_delim("../Blood_ATAC/Blood_ATAC/data/Waszak_2015/PU1_downloaded.txt", delim = " ", col_name = c("file_name"))

#Split into file names
file_names = dplyr::mutate(files, file_name = strsplit(URL, "/") %>% purrr::map_chr(., tail, 1) %>% unlist())
