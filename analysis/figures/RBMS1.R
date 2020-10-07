library("dplyr")
library("gwasvcf")
library("ggplot2")

#Visualise effect sizes
#' A general function to quickly import tabix indexed tab-separated files into data_frame
#'
#' @param tabix_file Path to tabix-indexed text file
#' @param param An instance of GRanges, RangedData, or RangesList
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to readr::read_delim()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export
scanTabixDataFrame <- function(tabix_file, param, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = readr::read_delim(result, delim = "\t", ...)[1,]
      }else{
        result = paste(x, collapse = "\n")
        result = readr::read_delim(result, delim = "\t", ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}

import_eQTLCatalogue <- function(ftp_path, region, selected_gene_id, column_names, verbose = TRUE){
  
  if(verbose){
    print(ftp_path)
  }
  
  #Fetch summary statistics with Rsamtools
  summary_stats = scanTabixDataFrame(ftp_path, region, col_names = column_names)[[1]] %>%
    dplyr::filter(gene_id == selected_gene_id)
  
  #Remove rsid duplicates and multi-allelic variant
  summary_stats = dplyr::select(summary_stats, -rsid) %>% 
    dplyr::distinct() %>% #rsid duplicates
    dplyr::mutate(id = paste(chromosome, position, sep = ":")) %>% 
    dplyr::group_by(id) %>% 
    dplyr::mutate(row_count = n()) %>% dplyr::ungroup() %>% 
    dplyr::filter(row_count == 1) #Multialllics
  
  return(summary_stats)
}



#Import summary stats
eqtl_data = readr::read_tsv("results/figure_data/BLUEPRINT_PE.T-cell_ge.nominal.sorted.tsv.gz") %>%
  dplyr::filter(molecular_trait_id == "ENSG00000153250") %>%
  dplyr::mutate(LP = -log(pvalue,10)) %>%
  dplyr::select(position, LP) %>%
  dplyr::mutate(type = "RBMS1 eQTL")
gwas_data = gwasvcf::query_gwas("results/figure_data/LC_GWAS_subset.GRCh38.sorted.vcf.gz", chrompos = "2:160100000-160700000") %>% 
  gwasvcf::vcf_to_granges() %>% 
  as.data.frame() %>%
  dplyr::as_tibble() %>%
  dplyr::transmute(position = start, LP, type = "Lymphocyte count")
gwas_data_selected = dplyr::semi_join(gwas_data, eqtl_data, by = "position")

joint_data = dplyr::bind_rows(eqtl_data, gwas_data)

#Import eQTL credible sets
cs = readr::read_tsv("results/figure_data/BLUEPRINT_PE.T-cell_ge.purity_filtered.txt.gz") %>%
  dplyr::filter(phenotype_id == "ENSG00000153250")

#Flag the credible set in GWAS results
gwas_flagged = dplyr::mutate(joint_data, in_cs = ifelse(position %in% cs$pos, TRUE, FALSE))

#Make manhattan plots
ggplot(eqtl_data, aes(y = -log(pvalue), x = position)) + geom_point()
manhattan = ggplot(gwas_flagged, aes(y = LP, x = position, color = in_cs)) + geom_point() + facet_grid(type~.)
ggsave("results/figures/RBMS1_manhattan.pdf", plot = manhattan, width = 5, height = 4)

#Import all coloc results
files_list = list.files("results/figure_data/LC-ebi-a-GCST004627_rnaseq/")
paths = paste0("results/figure_data/LC-ebi-a-GCST004627_rnaseq/", files_list)
path_list = setNames(as.list(paths), files_list)
coloc_res = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "dataset")
rbms1_colocs = dplyr::filter(coloc_res, molecular_trait_id == "ENSG00000153250", variant == "chr2_160468964_A_T") %>%
  dplyr::mutate(qtl_subset = stringr::str_replace(qtl_subset, pattern = "BLUEPRINT_SE", replacement = "BLUEPRINT")) %>%
  dplyr::mutate(qtl_subset = stringr::str_replace(qtl_subset, pattern = "BLUEPRINT_PE", replacement = "BLUEPRINT"))


#Import all coloc results
files_list = list.files("results/figure_data/LC-ebi-a-GCST004627_GTEx_V8/")
paths = paste0("results/figure_data/LC-ebi-a-GCST004627_GTEx_V8/", files_list)
path_list = setNames(as.list(paths), files_list)
coloc_res = purrr::map_df(path_list, ~readr::read_tsv(.), .id = "dataset")
gtex_colocs = dplyr::filter(coloc_res, molecular_trait_id == "ENSG00000153250", variant == "chr2_160468964_A_T")

rbms1_colocs = dplyr::bind_rows(rbms1_colocs, gtex_colocs)

#Visualise
coloc_plot = ggplot(rbms1_colocs, aes(x = PP.H3.abf, y = PP.H4.abf, label = qtl_subset)) + geom_point() + geom_text()
ggsave("results/figures/RBMS1_coloc_plot.pdf", plot = coloc_plot, width = 4, height = 4)



region_granges = GenomicRanges::GRanges(
  seqnames = "2", 
  ranges = IRanges::IRanges(start = 160468964-1, end = 160468964), 
  strand = "*")
region_granges


tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::as_tibble() %>%
  dplyr::filter(quant_method == "ge") %>%
  dplyr::mutate(qtl_id = paste0(study, ".", qtl_group, "_ge"))

#Import effect sizes
ftp_path_list = setNames(as.list(tabix_paths$ftp_path), tabix_paths$qtl_id)

#Extract column names from first file
column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))

#Import summmary stats
summary_list = purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region_granges, selected_gene_id = "ENSG00000153250", column_names))
summary_df = purrr::map_df(summary_list, identity, .id = "qtl_id")



imported_tabix_paths = read.delim("https://raw.githubusercontent.com/eQTL-Catalogue/eQTL-Catalogue-resources/master/tabix/tabix_ftp_paths_imported.tsv", sep = "\t", header = TRUE, stringsAsFactors = FALSE) %>% 
  dplyr::as_tibble() %>%
  dplyr::mutate(qtl_id = paste0(study, ".", qtl_group, "_ge"))
#Import effect sizes
ftp_path_list = setNames(as.list(imported_tabix_paths$ftp_path), imported_tabix_paths$qtl_id)

#Extract column names from first file
column_names = colnames(readr::read_tsv(ftp_path_list[[1]], n_max = 1))

#Import summmary stats
summary_list = purrr::map(ftp_path_list, ~import_eQTLCatalogue(., region_granges, selected_gene_id = "ENSG00000153250", column_names))
summary_df2 = purrr::map_df(summary_list, identity, .id = "qtl_id")
all_effects = dplyr::bind_rows(summary_df, summary_df2)



#Make an effect size plot
forest_plot = ggplot(all_effects, aes(x = qtl_id, y = beta, ymin = beta - se, ymax = beta + se)) + 
  geom_point() + 
  geom_errorbar(width = 0.1) + 
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_hline(yintercept = 0)
ggsave("results/figures/RBMS1_forest_plot.pdf", plot = forest_plot, width = 10, height = 5)

#Make an effect size plot
ggplot(rbms1_colocs, aes(x = qtl_subset, y = PP.H4.abf)) + 
  geom_bar(stat = "identity") +
  theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




