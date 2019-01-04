library("oligo")
library("dplyr")

#Import and normalise CEL files
cel_files = list.celfiles("array/cel_files/", listGzipped = TRUE) %>% file.path("array/cel_files/", .)
affyRaw = read.celfiles(cel_files)
eset <- rma(affyRaw)
saveRDS(eset, "ImmVar_hugene_eset.rds")
