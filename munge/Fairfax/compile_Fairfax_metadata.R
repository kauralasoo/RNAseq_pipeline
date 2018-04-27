library(tidyr)
library("dplyr")
library(readr)

pair1 = readr::read_csv("metadata/Fairfax/1_pair.csv") %>%
  dplyr::filter(!is.na(sample))
pair19 = readr::read_csv("metadata/Fairfax/19_pair.csv") %>%
  dplyr::transmute(File_name, sample = sample2, treatment) %>%
  dplyr::filter(!is.na(sample))
pair21 = readr::read_csv("metadata/Fairfax/21_pair.csv") %>%
  dplyr::transmute(File_name, sample = sample, treatment) %>%
  dplyr::filter(!is.na(sample))

merged_data = dplyr::bind_rows(pair1, pair19, pair21) %>%
  dplyr::rename(file_name = File_name) %>%
  dplyr::mutate(sample_id = paste0("S", sample, "_", treatment)) %>%
  dplyr::select(sample_id, everything())
