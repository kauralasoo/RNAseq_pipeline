meta = read.table("~/projects/SampleArcheology/studies/GEUVADIS/1000_genomes_sample_metadata.tsv", sep = "\t", header = T) %>%
  dplyr::as_tibble() %>%
  dplyr::filter(Superpopulation.code == "EUR", Sex == "female")
write.table(meta, "~/projects/SampleArcheology/studies/GEUVADIS/1000_genomes_EUR_female.tsv", sep = "\t", row.names = F, quote =F)


#Compare allele frequencies
cedar = eQTLUtils::importVariantInformation("~/projects/RNAseq_pipeline/results/X_impute/chrX.dose.filtered.variant_information.txt.gz") %>%
  dplyr::mutate(cedar_AF = AC/AN) %>%
  dplyr::select(pos, ref, alt, cedar_AF)
reference = eQTLUtils::importVariantInformation("~/projects/RNAseq_pipeline/results/X_impute/chrX.female.variant_information.txt.gz") %>%
  dplyr::mutate(ref_AF = AC/AN) %>%
  dplyr::select(pos, ref, alt, ref_AF)

comparison = dplyr::left_join(cedar, reference, by = c("pos","ref", "alt")) %>%
  dplyr::filter(!is.na(ref_AF))

smoothScatter(comparison$cedar_AF, comparison$ref_AF)
ggplot(comparison, aes(x = cedar_AF, y = ref_AF)) + geom_point()

cor(comparison$cedar_AF, comparison$ref_AF)
