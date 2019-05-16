library("dplyr")
library("knitr")
library("data.table")


convertAccessionToLink <- function(accession){
  if(accession %like% "^E-MTAB-"){
    string = paste0("[",accession,"](https://www.ebi.ac.uk/arrayexpress/experiments/",accession,")")
    return(string)
  }
  if(accession %like% "^E-GEUV-"){
    string = paste0("[",accession,"](https://www.ebi.ac.uk/arrayexpress/experiments/",accession,")")
    return(string)
  }
  if(accession %like% "^PRJEB"){
    string = paste0("[",accession,"](https://www.ebi.ac.uk/ena/data/view/",accession,")")
    return(string)
  }
  if(accession %like% "^EGAD"){
    string = paste0("[",accession,"](https://www.ebi.ac.uk/ega/datasets/",accession,")")
    return(string)
  }
  if(accession %like% "^GSE"){
    string = paste0("[",accession,"](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=",accession,")")
    return(string)
  }
  if(accession %like% "^phs"){
    string = paste0("[",accession,"](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=",accession,")")
    return(string)
  }
  return(accession)
}

string2string <- function(acc_string){
  acc_list = as.list(unlist(strsplit(acc_string, "; ")))
  acc_urls = purrr::map(acc_list, convertAccessionToLink)
  return(paste(acc_urls, collapse = "; "))
}

#Import data
table = read.table("metadata/website_table/List of datasets - Sheet1.csv", sep = ",", header = TRUE, stringsAsFactors = F) %>% 
  dplyr::as_tibble() %>% 
  dplyr::filter(status == "processed") %>%
  dplyr::filter(!(study_name %in% c("Raj_2014", "Ye_2018")))

#Markdown table
final_result = dplyr::mutate(table, study_text = paste0("[",study_name,"](",publication_doi,")"),
                             tissue, gxa_link = ifelse(gxa_accession != "", 
                                                       paste0("[",gxa_accession,"](https://www.ebi.ac.uk/gxa/experiments/", gxa_accession,")"), 
                                                       "")) %>%
  dplyr::group_by(study_name) %>%
  dplyr::mutate(rna_text = string2string(rna_accession), genotype_text = string2string(genotype_accession)) %>%
  dplyr::ungroup() %>%
  dplyr::select(study_text, tissue, stimulation, n_samples, n_donors, gxa_link, rna_text, genotype_text, experiment_type)


microarray_results = dplyr::filter(final_result, experiment_type == "microarray") %>% dplyr::select(-experiment_type)
write.table(microarray_results, "metadata/website_table/array_markdown_table.md", sep = " | ", row.names = F, col.names = F, quote = F)   

rnaseq_results = dplyr::filter(final_result, experiment_type == "RNA-seq") %>% dplyr::select(-experiment_type)
write.table(rnaseq_results, "metadata/website_table/rnaseq_markdown_table.md", sep = " | ", row.names = F, col.names = F, quote = F)   



