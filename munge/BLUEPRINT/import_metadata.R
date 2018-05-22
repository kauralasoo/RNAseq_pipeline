library("tidyr")
library("readr")
library("data.table")

#### Monocytes ####
#Sample files
sample_files = readr::read_tsv("metadata/BLUEPRINT/monocyte_sample_file.map", col_names = c("sample_id","id1", "file_name", "ega_file")) %>%
  dplyr::select(-id1) %>%
  dplyr::filter(!(file_name %like% ".bam.")) %>%
  tidyr::separate(file_name, c("file_name", "suffix"), ".cip") %>%
  dplyr::select(-suffix) %>%
  dplyr::filter(!(file_name %like% "mcgill"))

#Import file names
file_names = readr::read_tsv("metadata/BLUEPRINT/monocyte_fastq_files.txt", col_names = "fastq") %>%
  tidyr::separate(fastq, c("prefix","ega_id", "file_name"), "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-prefix)

#Export sample metadata
mono_sample_names = dplyr::left_join(sample_files, file_names, by = "file_name") %>% dplyr::select(sample_id, ega_id, fastq, ega_file) %>%
  dplyr::filter(!is.na(fastq))
write.table(mono_sample_names, "metadata/BLUEPRINT/monocyte_compiled_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Export snakemake fastq paths
snakemake = dplyr::transmute(mono_sample_names, sample_id, fastq) %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/BLUEPRINT/fastq/monocytes/", fastq,"]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakemake, "metadata/BLUEPRINT/monocytes_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



#Identify missing fastq files from EGA
all_files = read.table("metadata/BLUEPRINT/monocyte_EGA_list.txt", stringsAsFactors = FALSE)
setdiff(all_files$V3, sample_names$ega_file)


#### Neutrophils ####
#Sample files
sample_files = readr::read_tsv("metadata/BLUEPRINT/neutrophil_sample_file.map", col_names = c("sample_id","id1", "file_name", "ega_file")) %>%
  dplyr::select(-id1) %>%
  dplyr::filter(!(file_name %like% ".bam.")) %>%
  tidyr::separate(file_name, c("file_name", "suffix"), ".cip") %>%
  dplyr::select(-suffix) %>%
  dplyr::filter(!(file_name %like% "mcgill")) %>%
  dplyr::filter(!(file_name %like% "McGill"))

#Import file names
file_names = readr::read_tsv("metadata/BLUEPRINT/neutrophil_fastq_files.txt", col_names = "fastq") %>%
  tidyr::separate(fastq, c("prefix","ega_id", "file_name"), "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-prefix)

neutro_sample_names = dplyr::left_join(sample_files, file_names, by = "file_name") %>% dplyr::select(sample_id, ega_id, fastq, ega_file) %>%
  dplyr::filter(!is.na(fastq))
write.table(neutro_sample_names, "metadata/BLUEPRINT/neutrophil_compiled_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)

#Export snakemake fastq paths
snakemake = dplyr::transmute(neutro_sample_names, sample_id, fastq) %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/BLUEPRINT/fastq/neutrophils/", fastq,"]")) %>% 
  dplyr::select(snakemake_string)
write.table(snakemake, "metadata/BLUEPRINT/neutrophils_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)


#Identify missing fastq files from EGA
all_files = read.table("metadata/BLUEPRINT/neutrophil_EGA_list.txt", stringsAsFactors = FALSE)
setdiff(all_files$V3, neutro_sample_names$ega_file)



#T-cell data
#Sample files
sample_files = readr::read_tsv("metadata/BLUEPRINT/Tcell_sample_file.map", col_names = c("sample_id","id1", "file_name", "ega_file")) %>%
  dplyr::select(-id1) %>%
  dplyr::filter(!(file_name %like% ".bam.")) %>%
  tidyr::separate(file_name, c("file_name", "suffix"), ".cip") %>%
  dplyr::select(-suffix) %>%
  dplyr::mutate(file_name = str_replace_all(file_name, "/", "_")) %>%
  dplyr::mutate(file_name = str_replace(file_name, "_McGill", "McGill"))


#Import file names
file_names = readr::read_tsv("metadata/BLUEPRINT/Tcell_fastq_files.txt", col_names = "fastq") %>%
  tidyr::separate(fastq, c("prefix","ega_id", "file_name"), "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-prefix)

tcell_sample_names = dplyr::left_join(sample_files, file_names, by = "file_name") %>% 
  dplyr::select(sample_id, ega_id, fastq, ega_file) %>%
  dplyr::filter(!is.na(fastq))





#Identify missing fastq files from EGA
all_files = read.table("metadata/BLUEPRINT/Tcell_EGA_list.txt", stringsAsFactors = FALSE)
setdiff(all_files$V3, neutro_sample_names$ega_file)

