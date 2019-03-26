library("dplyr")
library("purrr")
library("tidyr")
library("devtools")
library("optparse")
library("GenomicRanges")
load_all("../eQTLUtils/")

option_list <- list(
  #TODO look around if there is a package recognizing delimiter in dataset
  optparse::make_option(c("-p", "--pairs"), type="character", default="results/lead_vars.txt",
                        help="Pairs of phenotype_ids and snp_ids to be extracted from summary files. Tab separated", metavar = "type"),
  optparse::make_option(c("-q", "--qtlgroups"), type="character", default="results/qtl_groups.txt",
                        help="Sample metadata file path of genes used in expression-matrix. Tab separated", metavar = "type"),
  optparse::make_option(c("-o", "--outfile"), type="character", default="./extracted_variants.txt",
                        help="Path to the output file.", metavar = "type")
  )

opt <- optparse::parse_args(optparse::OptionParser(option_list=option_list))

#Import data
selected_pairs = read.table(opt$p, sep = "\t", stringsAsFactors = F, header = T)
qtl_groups = read.table(opt$q, sep = "\t", stringsAsFactors = F, header = F)[,1]
tabix_paths = setNames(paste0(qtl_groups, ".nominal.sorted.txt.gz"), qtl_groups)

#Fetch all pairs in all datasets
datasets = purrr::map_df(tabix_paths, ~eQTLUtils::qtltoolsTabixFetchPhenotypeVariantPairs(selected_pairs, .), .id = "qtl_group") %>%
  dplyr::mutate(z_score = abs(qnorm(p_nominal/2, mean = 0, sd = 1))*sign(beta)) %>%
  dplyr::mutate(se = beta/abs(z_score))

#Export results
write.table(datasets, opt$o, quote = F, row.names = F, sep = "\t")


