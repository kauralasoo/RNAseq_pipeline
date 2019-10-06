monocyte = read.table("results/mixupmapper/monocyte_CD14/MixupMapper/BestMatchPerTrait.txt", 
                      header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
bcell = read.table("results/mixupmapper/B-cell_CD19/MixupMapper/BestMatchPerTrait.txt", 
                      header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
neutrophil = read.table("results/mixupmapper/neutrophil_CD15/MixupMapper/BestMatchPerTrait.txt", 
                      header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
tcd4 = read.table("results/mixupmapper/T-cell_CD8/MixupMapper/BestMatchPerTrait.txt", 
                      header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
tcd8 = read.table("results/mixupmapper/T-cell_CD4/MixupMapper/BestMatchPerTrait.txt", 
                      header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
platelet = read.table("results/mixupmapper/platelet/MixupMapper/BestMatchPerTrait.txt", 
                  header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
rectum = read.table("results/mixupmapper/rectum/MixupMapper/BestMatchPerTrait.txt", 
                  header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
transverse = read.table("results/mixupmapper/transverse_colon/MixupMapper/BestMatchPerTrait.txt", 
                  header =T, stringsAsFactors = F) %>% dplyr::as_tibble()
ileum = read.table("results/mixupmapper/ileum/MixupMapper/BestMatchPerTrait.txt", 
                        header =T, stringsAsFactors = F) %>% dplyr::as_tibble()

all_results = dplyr::bind_rows(monocyte, bcell, neutrophil, tcd4, tcd8, platelet, rectum, transverse, ileum) %>%
  dplyr::filter(Mixup == "true") %>%
  dplyr::mutate(score_diff = abs(BestMatchingGenotypeScore) - abs(OriginalLinkedGenotypeScore)) %>%
  dplyr::arrange(score_diff) %>%
  dplyr::filter(score_diff > 0.5)

#Import current metadata 
current_meta = read.table("../SampleArcheology/studies/cleaned/CEDAR.tsv", sep = "\t", 
                          header = T, stringsAsFactors = F) %>%
  dplyr::as_tibble()
sex_df = dplyr::select(current_meta, genotype_id, sex) %>% distinct()

fixed_genotypes = dplyr::transmute(all_results, sample_id = Trait, 
                 genotype_id = OriginalLinkedGenotype, new_genotype_id = BestMatchingGenotype) %>%
  dplyr::left_join(sex_df, by = c("new_genotype_id" = "genotype_id")) %>%
  dplyr::transmute(sample_id, new_genotype_id, new_sex = sex)

new_meta = dplyr::left_join(current_meta, fixed_genotypes, by = "sample_id") %>%
  dplyr::mutate(genotype_id = ifelse(!is.na(new_genotype_id), new_genotype_id, genotype_id)) %>%
  dplyr::mutate(sex = ifelse(!is.na(new_sex), new_sex, sex)) %>%
  dplyr::select(-new_genotype_id, -new_sex)
write.table(new_meta, "../SampleArcheology/studies/cleaned/CEDAR.tsv", sep = "\t", row.names = F, quote = F)





