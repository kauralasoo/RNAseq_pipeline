library("tidyr")
library("readr")
library("data.table")

#Sample files
sample_files = readr::read_tsv("metadata/BLUEPRINT/monocyte_sample_file.map", col_names = c("sample_id","id1", "file_name", "id2")) %>%
  dplyr::select(-id1, -id2) %>%
  dplyr::filter(!(file_name %like% ".bam.")) %>%
  tidyr::separate(file_name, c("file_name", "suffix"), ".cip") %>%
  dplyr::select(-suffix) %>%
  dplyr::filter(!(file_name %like% "mcgill"))

#Import file names
file_names = readr::read_tsv("metadata/BLUEPRINT/monocyte_fastq_files.txt", col_names = "fastq") %>%
  tidyr::separate(fastq, c("prefix","ega_id", "file_name"), "_", extra = "merge", remove = FALSE) %>%
  dplyr::select(-prefix)

sample_names = dplyr::left_join(sample_files, file_names, by = "file_name") %>% dplyr::select(sample_id, ega_id, fastq) %>%
  dplyr::filter(!is.na(fastq))
write.table(sample_names, "metadata/BLUEPRINT/monocyte_compiled_metadata.txt", sep = "\t", row.names = FALSE, quote = FALSE)
