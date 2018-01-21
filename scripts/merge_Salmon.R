library("readr")
library("optparse")
library("tximport")

#Parse command-line options
option_list <- list(
  make_option(c("-s", "--samples"), type="character", default=NULL,
              help="Comma-separated list of sample ids.", metavar = "type"),
  make_option(c("-d", "--dir"), type="character", default=NULL,
              help="Directory of the featureCounts output files.", metavar = "type"),
  make_option(c("-o", "--out"), type="character", default=NULL,
              help="Path to the output file.", metavar = "type")
)
opt <- parse_args(OptionParser(option_list=option_list))

#Test input
#opt = list(s = "ERR188067,ERR188081,ERR188084", d = "processed/GEUVADIS/salmon/Ensembl_87", o = "processed/GEUVADIS/matrices/Ensembl_87.salmon_txrevise.rds")

#Extract sample ids
sample_ids = unlist(strsplit(opt$s, ","))

#Construct file names
file_names = setNames(file.path(opt$d, sample_ids, "quant.sf"), sample_ids)

#Import quant results
tx_abundances = tximport(file_names, type = "salmon", txOut = TRUE, importer = read_tsv, dropInfReps = TRUE, ignoreTxVersion = FALSE)
saveRDS(tx_abundances, opt$o)


