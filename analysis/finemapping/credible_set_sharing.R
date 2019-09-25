library("ggplot2")
library("tidyr")

importFinemapResults <- function(finemap_path){
  
  short_paths = list.files(finemap_path)
  full_paths = list.files(finemap_path, full.names = T)
  
  #Extract qtl_group
  paths_df = dplyr::data_frame(file_name = full_paths, short_name = short_paths) %>%
    tidyr::separate(short_name, c("study","qtl_group", "quant","batch", "suffix"), sep = "\\.") %>%
    dplyr::select(study, qtl_group, quant, file_name) %>%
    dplyr::distinct() 
  
  all_files = setNames(full_paths, full_paths)
  imported_files = purrr::map_df(all_files, ~read.table(., header = T, stringsAsFactors = F) %>% dplyr::as_tibble(), .id = "file_name") %>%
    dplyr::left_join(paths_df, by = c("file_name")) %>%
    dplyr::select(-file_name)
}

calculateCSSharing <- function(cs_list1, cs_list2){
  cs1_df = dplyr::mutate(cs_list1, full_cs_id = paste(phenotype_id, cs_id, sep = "_")) %>%
    dplyr::select(full_cs_id, variant_id) %>%
    dplyr::group_by(full_cs_id) %>%
    dplyr::mutate(cs_size1 = n()) %>%
    dplyr::filter(cs_size1 < 100) %>%
    dplyr::ungroup() %>%
    dplyr::rename(full_cs_id1 = full_cs_id)
  cs2_df = dplyr::mutate(cs_list2, full_cs_id = paste(phenotype_id, cs_id, sep = "_"))  %>%
    dplyr::select(full_cs_id, variant_id) %>%
    dplyr::group_by(full_cs_id) %>%
    dplyr::mutate(cs_size2 = n()) %>%
    dplyr::filter(cs_size2 < 100) %>%
    dplyr::ungroup() %>%
    dplyr::rename(full_cs_id2 = full_cs_id)
  
  #Join the two datsets
  overlaps = dplyr::left_join(cs1_df, cs2_df, by = "variant_id") %>%
    dplyr::group_by(full_cs_id1, full_cs_id2) %>%
    dplyr::mutate(overlap_size = n()) %>%
    dplyr::mutate(overlap_fraction = overlap_size/min(cs_size1, cs_size2)) %>%
    dplyr::filter(overlap_fraction > 0.2) %>%
    dplyr::ungroup()
  
  trait1_rep = length(unique(overlaps$full_cs_id1))/length(unique(cs1_df$full_cs_id1))
  trait2_rep = length(unique(overlaps$full_cs_id2))/length(unique(cs2_df$full_cs_id2))
  result = dplyr::tibble(sharing1 = trait1_rep, sharing2 = trait2_rep)
  return(result)
}

#Import finemapping results for different studies
geuvadis_results = importFinemapResults("results/finemap/GEUVADIS/txt/") %>% 
  dplyr::filter(quant == "gene_counts")

gencord_results = importFinemapResults("results/finemap/GENCORD/txt/") %>% 
  dplyr::filter(quant == "gene_counts")

twinsuk_results = importFinemapResults("results/finemap/TwinsUK/txt/") %>% 
  dplyr::filter(quant == "gene_counts")

lepik_2017_results = importFinemapResults("results/finemap/Lepik_2017/txt/") %>% 
  dplyr::filter(quant == "gene_counts")

#Make list
full_df = dplyr::bind_rows(twinsuk_results, gencord_results, geuvadis_results, lepik_2017_results) %>% 
  dplyr::group_by(study, qtl_group)
keys = dplyr::group_keys(full_df) %>%
  dplyr::mutate(list_label = paste(qtl_group, study, sep = "_"))
list = dplyr::group_split(full_df)
names(list) = keys$list_label

#Find all pairs of traits
pairs = t(combn(keys$list_label, 2)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  dplyr::as_tibble() %>%
  dplyr::rename(trait1 = V1, trait2 = V2)

#Caclulate sharing
list1 = list[pairs$trait1]
list2 = list[pairs$trait2]
sharing_list = purrr::map2(list1, list2, ~calculateCSSharing(.x, .y))
sharing_df = purrr::map_df(sharing_list, identity) %>%
  dplyr::bind_cols(pairs, .) %>%
  dplyr::mutate(mean_sharing = (sharing1+sharing2)/2)

sharing_result = dplyr::bind_rows(dplyr::transmute(sharing_df, discovery = trait1, replication = trait2, replication_rate = mean_sharing),
dplyr::transmute(sharing_df, discovery = trait2, replication = trait1, replication_rate = mean_sharing), 
dplyr::transmute(keys, discovery = list_label, replication = list_label, replication_rate = 1))

ggplot(sharing_result, aes(y = discovery, x = replication, fill = replication_rate)) + 
  geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", 
                       mid = "#FFFFBF", high = "#E24C36", name = "Relative effect", midpoint = 0.5, limits = c(0,1))


sharing_df2 = tidyr::pivot_wider(sharing_result, id_cols = "discovery", names_from = "replication", values_from = "replication_rate")
sharing_mat = dplyr::select(sharing_df2, -discovery) %>% as.matrix()
rownames(sharing_mat) = sharing_df2$discovery
pdf("results/figures/rnaseq_qtl_sharing.pdf", width = 7, height = 7)
pheatmap(sharing_mat)
dev.off()





#Redo the same analysis for array datasets
#Import finemapping results for different studies
kasela_2017 = importFinemapResults("results/finemap/Kasela_2017/txt/")

cedar = importFinemapResults("results/finemap/CEDAR/txt/") 

fairfax_2014 = importFinemapResults("results/finemap/Fairfax_2014//txt/") 

fairfax_2012 = importFinemapResults("results/finemap/Fairfax_2012//txt/") 

naranbhai_2015 = importFinemapResults("results/finemap/Naranbhai_2015//txt/")


#Make list
full_df = dplyr::bind_rows(kasela_2017, cedar, fairfax_2014, fairfax_2012, naranbhai_2015) %>% 
  dplyr::group_by(study, qtl_group)
keys = dplyr::group_keys(full_df) %>%
  dplyr::mutate(list_label = paste(qtl_group, study, sep = "_"))
list = dplyr::group_split(full_df)
names(list) = keys$list_label

#Find all pairs of traits
pairs = t(combn(keys$list_label, 2)) %>% 
  as.data.frame(stringsAsFactors = FALSE) %>% 
  dplyr::as_tibble() %>%
  dplyr::rename(trait1 = V1, trait2 = V2)

#Caclulate sharing
list1 = list[pairs$trait1]
list2 = list[pairs$trait2]
sharing_list = purrr::map2(list1, list2, ~calculateCSSharing(.x, .y))
sharing_df = purrr::map_df(sharing_list, identity) %>%
  dplyr::bind_cols(pairs, .) %>%
  dplyr::mutate(mean_sharing = (sharing1+sharing2)/2)

sharing_result = dplyr::bind_rows(dplyr::transmute(sharing_df, discovery = trait1, replication = trait2, replication_rate = mean_sharing),
                                  dplyr::transmute(sharing_df, discovery = trait2, replication = trait1, replication_rate = mean_sharing), 
                                  dplyr::transmute(keys, discovery = list_label, replication = list_label, replication_rate = 1))

ggplot(sharing_result, aes(y = discovery, x = replication, fill = replication_rate)) + 
  geom_tile() + 
  scale_fill_gradient2(space = "Lab", low = "#4575B4", 
                       mid = "#FFFFBF", high = "#E24C36", name = "CS sharing", midpoint = 0.5, limits = c(0,1))


sharing_df2 = tidyr::pivot_wider(sharing_result, id_cols = "discovery", names_from = "replication", values_from = "replication_rate")
sharing_mat = dplyr::select(sharing_df2, -discovery) %>% as.matrix()
rownames(sharing_mat) = sharing_df2$discovery
pdf("results/figures/array_qtl_sharing.pdf", width = 7, height = 7)
pheatmap(sharing_mat)
dev.off()


