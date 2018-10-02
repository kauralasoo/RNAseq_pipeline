library("SeqArray")
library("SNPRelate")

#Convert to GDS
snpgdsVCF2GDS("vcf/TwinsUK_GRCh38.annotated.vcf.gz", "vcf/TwinsUK_GRCh38.annotated.gds", method="biallelic.only")


#Perform LD pruning
genofile <- snpgdsOpen("vcf/TwinsUK_GRCh38.annotated.gds")

# Try different LD thresholds for sensitivity analysis
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, ld.threshold=0.2)

# Get all selected snp id
snpset.id <- unlist(unname(snpset))
head(snpset.id)

#Calculate IBD
ibd <- snpgdsIBDMoM(genofile, snp.id=snpset.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

#Identify twins
mz_twins = dplyr::filter(ibd.coeff, k0 < 0.2, k1 < 0.2) %>%
  dplyr::as_tibble() %>%
  dplyr::select(ID1, ID2, kinship) %>%
  dplyr::mutate(twin_id = paste0("MZ_", 1:length(ID2))) %>%
  tidyr::gather("twin_name", "genotype_id", ID1:ID2) %>%
  dplyr::arrange(twin_id)
dz_twins = dplyr::filter(ibd.coeff, k0 > 0.1, k1 > 0.25) %>%
  dplyr::as_tibble() %>%
  dplyr::select(ID1, ID2, kinship) %>%
  dplyr::mutate(twin_id = paste0("DZ_", 1:length(ID2))) %>%
  tidyr::gather("twin_name", "genotype_id", ID1:ID2) %>%
  dplyr::arrange(twin_id)
black_list = dplyr::filter(ibd.coeff, k1 < 0.2, k0 > 0.5, k0 < 0.8) %>%  
  dplyr::select(ID1, ID2, kinship) %>%
  tidyr::gather("twin_name", "genotype_id", ID1:ID2)

#Merge together
twins_df = dplyr::bind_rows(mz_twins, dz_twins) %>%
  dplyr::select(genotype_id, everything())

#Identifty families manually
dplyr::group_by(dz_twins, genotype_id) %>% 
  dplyr::mutate(count = length(genotype_id)) %>% 
  dplyr::filter(count > 1) %>% 
  dplyr::arrange(genotype_id) %>% View()
fam1 = c("DZ_50", "DZ_51", "DZ_52", "DZ_53", "DZ_54", "DZ_142")
fam2 = c("DZ_158", "DZ_159", "DZ_160", "DZ_172", "DZ_170")

dz_unrelated = dplyr::mutate(dz_twins, twin_id = ifelse(twin_id %in% fam1, "DZ_FAM_1", twin_id)) %>%
  dplyr::mutate(twin_id = ifelse(twin_id %in% fam2, "DZ_FAM_1", twin_id)) %>%
  dplyr::select(genotype_id, twin_id) %>%
  dplyr::group_by(twin_id) %>%
  dplyr::arrange(genotype_id) %>%
  dplyr::filter(row_number() == 1) %>%
  dplyr::ungroup()

#Unrelated 
all_individuals = dplyr::as_tibble(ibd.coeff) %>% 
  dplyr::select(ID1) %>%
  dplyr::ungroup() %>%
  dplyr::distinct() %>%
  dplyr::rename(genotype_id = ID1)
no_twins = dplyr::anti_join(all_individuals, twins_df) %>%
  dplyr::anti_join(black_list) %>%
  dplyr::distinct()

#Identify all unrelated individuals
unrelated = c(no_twins$genotype_id, dz_unrelated$genotype_id, dplyr::filter(mz_twins, twin_name == "ID1")$genotype_id)
write.table(unrelated, "metadata/TwinsUK/TwinsUK_unrelated_individuals.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

#Perform LD pruning
genofile <- snpgdsOpen("vcf/TwinsUK_GRCh38.annotated.gds")

#Extract samples
sample_ids = read.gdsn(index.gdsn(genofile, "sample.id"))

# Try different LD thresholds for sensitivity analysis
set.seed(1000)
snpset <- snpgdsLDpruning(genofile, sample.id = unrelated, ld.threshold=0.2)

#Calculate IBD
ibd <- snpgdsIBDMoM(genofile, sample.id = unrelated, snp.id=snpset.id,
                    maf=0.05, missing.rate=0.05, num.thread=2)
ibd.coeff <- snpgdsIBDSelection(ibd)
head(ibd.coeff)

plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="YRI samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)


# Perform PCA
pca <- snpgdsPCA(genofile, snp.id = snpset.id, sample.id = unrelated, num.thread=2)

tab <- data.frame(sample.id = pca$sample.id,
                  EV1 = pca$eigenvect[,1],    # the first eigenvector
                  EV2 = pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(tab)
plot(tab$EV2, tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1")

