library("SNPRelate")
library("GDSArray")
library("devtools")
load_all("../eQTLUtils/")
load_all("../seqUtils/")
library("ggplot2")
library("cqn")

#Make a GDS file
snpgdsVCF2GDS("results/genotypes/GENCORD/GENCORD_GRCh38.filtered.vcf.gz", 
              "results/genotypes/GENCORD/GENCORD_GRCh38.filtered.gds", 
              method = "copy.num.of.ref")
snpgdsVCF2GDS("results/genotypes/Fairfax_2014/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz", 
              "results/genotypes/Fairfax_2014/Fairfax_2014_GRCh38.filtered.renamed.gds", 
              method = "copy.num.of.ref")
snpgdsVCF2GDS("results/genotypes/CEDAR/CEDAR_GRCh38.filtered.renamed.vcf.gz", 
              "results/genotypes/CEDAR/CEDAR_GRCh38.filtered.renamed.gds", 
              method = "copy.num.of.ref")
snpgdsVCF2GDS("results/genotypes/TwinsUK/TwinsUK_GRCh38.filtered.vcf.gz", 
              "results/genotypes/TwinsUK/TwinsUK_GRCh38.filtered.gds", 
              method = "copy.num.of.ref")

#Import GENCORD lead variants
GENCORD_LCL = eQTLUtils::importQTLtoolsTable("processed/GENCORD/qtltools/output/featureCounts/LCL.permuted.txt.gz")

#import variant info
gds_file = "results/genotypes/GENCORD/GENCORD_GRCh38.filtered.gds"
var_info = importVariantInformationFromGDS(gds_file)

#Selected QTL
qtl = dplyr::filter(GENCORD_LCL, phenotype_id == "ENSG00000166750")
selected_snp_id = qtl$snp_id

#Extract genotype
gencord_se = readRDS("results/SummarizedExperiments/GENCORD.rds") %>%
  filterSummarizedExperiment(filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
cqn_se = normaliseSE_cqn(gencord_se)

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, cqn_se, "cqn")
selected_genotype = extractVariantGenotypeFromGDS(selected_snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = cell_type, gene_name = "SLFN5")

qtl_plot = plotQtlRow(selected_data)
ggsave("results/figures/SLFN5_GENCORD.pdf", width = 4, height = 4)

#In TwinsUK
#import variant info
gds_file = "results/genotypes/TwinsUK/TwinsUK_GRCh38.filtered.gds"
var_info = importVariantInformationFromGDS(gds_file)

#Extract genotype
twinsuk_se = readRDS("results/SummarizedExperiments/TwinsUK.rds") %>%
  filterSummarizedExperiment(filter_rna_qc = TRUE, filter_genotype_qc = TRUE)
cqn_se = normaliseSE_cqn(twinsuk_se)

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, cqn_se, "cqn")
selected_genotype = extractVariantGenotypeFromGDS(selected_snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = cell_type, gene_name = "SLFN5")

qtl_plot = plotQtlRow(selected_data)
ggsave("results/figures/SLFN5_TwinsUK.pdf", width = 5, height = 4)


#Import fairfax lead variants
FF_LPS = eQTLUtils::importQTLtoolsTable("processed/Fairfax_2014/qtltools/output/array/monocyte_LPS2.permuted.txt.gz") 

#Selected QTL
qtl = dplyr::filter(FF_LPS, group_id == "ENSG00000166750")
fairfax_se = readRDS("results/SummarizedExperiments/array_norm/Fairfax_2014.rds")
assays(fairfax_se) = assays(fairfax_se)[c(1,3)]

#import variant info
gds_file = "results/genotypes/Fairfax_2014/Fairfax_2014_GRCh38.filtered.renamed.gds"
var_info = importVariantInformationFromGDS(gds_file)

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, fairfax_se, "norm_exprs")
selected_genotype = extractVariantGenotypeFromGDS(selected_snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = condition, gene_name = "SLFN5")
qtl_plot = plotQtlRow(selected_data)
ggsave("results/figures/SLFN5_Fairfax_2014.pdf",qtl_plot, width = 4, height = 4)


### CEDAR ###
cedar_se = readRDS("results/SummarizedExperiments/array_norm/CEDAR.rds")

gds_file = "results/genotypes/CEDAR/CEDAR_GRCh38.filtered.renamed.gds"
var_info = importVariantInformationFromGDS(gds_file)

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, cedar_se, "norm_exprs")
selected_genotype = extractVariantGenotypeFromGDS(selected_snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = qtl_group, gene_name = "SLFN5")
qtl_plot = plotQtlRow(selected_data)
plot = ggplot2::ggplot(selected_data, ggplot2::aes(x = genotype_text, y = norm_exp, color = condition_name)) + 
  ggplot2::facet_wrap(~condition_name, scales = "fixed", nrow = 1) + 
  ggplot2::geom_boxplot(outlier.shape = NA) + 
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
  ggplot2::ylab(paste0(selected_data$gene_name[1], " expression")) +
  ggplot2::xlab(selected_data$snp_id[1]) + 
  ggplot2::theme_light() + 
  theme(strip.text.x = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))


#MGAT3
CEDAR_B = eQTLUtils::importQTLtoolsTable("processed/CEDAR/qtltools/output/array/B-cell_CD19.permuted.txt.gz")
qtl = dplyr::filter(CEDAR_B, group_id == "ENSG00000128268")

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, cedar_se, "norm_exprs")
selected_genotype = extractVariantGenotypeFromGDS(qtl$snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = qtl_group, gene_name = "MGAT3")
qtl_plot = plotQtlRow(selected_data)
plot = ggplot2::ggplot(selected_data, ggplot2::aes(x = genotype_text, y = norm_exp, color = condition_name)) + 
  ggplot2::facet_wrap(~condition_name, scales = "fixed", nrow = 1) + 
  ggplot2::geom_boxplot(outlier.shape = NA) + 
  ggplot2::geom_jitter(position = ggplot2::position_jitter(width = .2), size = 0.5) + 
  ggplot2::ylab(paste0(selected_data$gene_name[1], " expression")) +
  ggplot2::xlab(selected_data$snp_id[1]) + 
  ggplot2::theme_light() + 
  theme(strip.text.x = element_text(colour = "grey10"), strip.background = element_rect(fill = "grey85"))
ggsave("results/figures/MGAT3_CEDAR.pdf", plot = plot, width = 10, height = 6)


#Replicate in Fairfax
fairfax_2012_se = readRDS("results/SummarizedExperiments/array_norm/Fairfax_2012.rds")
fairfax_all = cbind(fairfax_se, fairfax_2012_se)
fairfax_naive = fairfax_all[,fairfax_all$condition == "naive"]

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, fairfax_naive, "norm_exprs")
selected_genotype = extractVariantGenotypeFromGDS(qtl$snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = qtl_group, gene_name = "MGAT3")
qtl_plot = plotQtlRow(selected_data)
ggsave("results/figures/MGAT3_Fairfax_2012.pdf", plot = qtl_plot, width = 4, height = 4)

