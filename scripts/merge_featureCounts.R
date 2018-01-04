library("readr")
library("optparse")

#' Import featureCounts output table into R.
#'
#' Skips the first comment line.
#' 
#' @param sample_dir Path to the directory containing the count files.
#' @param sample_names Vector of sample names.
#' @param counts_suffix Suffix of the counts file.
#' @param sub_dir If TRUE, count files are nested in subfolders named after sample names, otherwise counts files
#' are directly in sample_dir.
#' @return Counts matrix where the first two columns are gene_id and feature length.
#' @author Kaur Alasoo
#' @export 
loadCounts <- function(sample_dir, sample_names, counts_suffix = ".counts.txt", sub_dir = TRUE){
  #Load featureCounts output into R
  matrix = c()
  for (i in c(1:length(sample_names))){
    if (sub_dir == TRUE){
      path = file.path(sample_dir, sample_names[i], paste(sample_names[i], counts_suffix, sep = ""))
    } else {
      path = file.path(sample_dir, paste(sample_names[i], counts_suffix, sep = ""))      
    }
    print(sample_names[i])
    table = readr::read_tsv(path, skip = 1, col_types = "cccccii")
    print(head(table))
    if (i == 1){
      matrix = table[,c(1,6,7)]
    }
    else{
      matrix = cbind(matrix, table[,7])
    }
  }
  colnames(matrix) = c("gene_id", "length", sample_names)
  return(matrix)
}

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
#opt = list(s = "ERR188067,ERR188081,ERR188284", d = "processed/GEUVADIS/featureCounts/", o = "processed/GEUVADIS/matrices/featureCounts.txt")

#Extract sample ids
sample_ids = unlist(strsplit(opt$s, ","))
count_matrix = loadCounts(opt$d, sample_ids, sub_dir = FALSE, counts_suffix = ".featureCounts.txt")
write.table(count_matrix, opt$o, row.names = FALSE, quote = FALSE, sep = "\t")


