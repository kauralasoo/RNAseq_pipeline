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




#### paired end data
#### Neutrophils ####
#Sample files
sample_files = readr::read_tsv("metadata/BLUEPRINT/neutrophil_sample_file.map", col_names = c("sample_id","id1", "file_name", "ega_file")) %>%
  dplyr::select(-id1) %>%
  dplyr::filter(!(file_name %like% ".bam.")) %>%
  tidyr::separate(file_name, c("file_name", "suffix"), ".cip") %>%
  dplyr::select(-suffix) %>%
  dplyr::filter((file_name %like% "mcgill") | (file_name %like% "McGill")) %>%
  dplyr::mutate(file_name = str_replace_all(file_name, "/", "_")) %>%
  dplyr::mutate(file_name = str_replace(file_name, "_McGill", "McGill"))

#Import file names
file_names = readr::read_tsv("metadata/BLUEPRINT/neutrophil_fastq_files.txt", col_names = "fastq") %>%
  tidyr::separate(fastq, c("prefix","ega_id", "file_name"), "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-prefix)

neutro_sample_names = dplyr::left_join(sample_files, file_names, by = "file_name") %>% dplyr::select(sample_id, ega_id, fastq, ega_file) %>%
  dplyr::filter(!is.na(fastq))

#Make snakemale string
fq1 = dplyr::transmute(neutro_sample_names, sample_id, fq1 = fastq) %>% dplyr::filter(fq1 %like% "pair1")
fq2 = dplyr::transmute(neutro_sample_names, sample_id, fq2 = fastq) %>% dplyr::filter(fq2 %like% "pair2")
neutrophil_snakemake = dplyr::left_join(fq1, fq2, by = "sample_id") %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/BLUEPRINT/fastq/neutrophils/", fq1, ", processed/BLUEPRINT/fastq/neutrophils/", fq2,"]")) %>% 
  dplyr::select(snakemake_string)

#### Neutrophils ####
#Sample files
sample_files = readr::read_tsv("metadata/BLUEPRINT/monocyte_sample_file.map", col_names = c("sample_id","id1", "file_name", "ega_file")) %>%
  dplyr::select(-id1) %>%
  dplyr::filter(!(file_name %like% ".bam.")) %>%
  tidyr::separate(file_name, c("file_name", "suffix"), ".cip") %>%
  dplyr::select(-suffix) %>%
  dplyr::filter((file_name %like% "mcgill") | (file_name %like% "McGill")) %>%
  dplyr::mutate(file_name = str_replace_all(file_name, "/", "_")) %>%
  dplyr::mutate(file_name = str_replace(file_name, "_McGill", "McGill"))

#Import file names
file_names = readr::read_tsv("metadata/BLUEPRINT/monocyte_fastq_files.txt", col_names = "fastq") %>%
  tidyr::separate(fastq, c("prefix","ega_id", "file_name"), "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-prefix)

neutro_sample_names = dplyr::left_join(sample_files, file_names, by = "file_name") %>% dplyr::select(sample_id, ega_id, fastq, ega_file) %>%
  dplyr::filter(!is.na(fastq))

#Make snakemale string
fq1 = dplyr::transmute(neutro_sample_names, sample_id, fq1 = fastq) %>% dplyr::filter(fq1 %like% "pair1")
fq2 = dplyr::transmute(neutro_sample_names, sample_id, fq2 = fastq) %>% dplyr::filter(fq2 %like% "pair2")
monocyte_snakemake = dplyr::left_join(fq1, fq2, by = "sample_id") %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/BLUEPRINT/fastq/monocytes/", fq1, ", processed/BLUEPRINT/fastq/monocytes/", fq2,"]")) %>% 
  dplyr::select(snakemake_string)


#### T-cell data ####
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

#Make snakemale string
fq1 = dplyr::transmute(tcell_sample_names, sample_id, fq1 = fastq) %>% dplyr::filter(fq1 %like% "pair1")
fq2 = dplyr::transmute(tcell_sample_names, sample_id, fq2 = fastq) %>% dplyr::filter(fq2 %like% "pair2")
tcell_snakemake = dplyr::left_join(fq1, fq2, by = "sample_id") %>% 
  dplyr::mutate(snakemake_string = paste0(sample_id,": [", "processed/BLUEPRINT/fastq/T-cells/", fq1, ", processed/BLUEPRINT/fastq/T-cells/", fq2,"]")) %>% 
  dplyr::select(snakemake_string)

all_snakemake = dplyr::bind_rows(monocyte_snakemake, neutrophil_snakemake, tcell_snakemake)
write.table(all_snakemake, "metadata/BLUEPRINT/PE_snakemake_string.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



#Identify missing fastq files from EGA
all_files = read.table("metadata/BLUEPRINT/Tcell_EGA_list.txt", stringsAsFactors = FALSE)
setdiff(all_files$V3, neutro_sample_names$ega_file)

