library("SNPRelate")
library("GDSArray")
library("devtools")
load_all("../eQTLUtils/")
load_all("../seqUtils/")
library("ggplot2")
library("cqn")
library("SummarizedExperiment")

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
FF_naive = eQTLUtils::importQTLtoolsTable("processed/Fairfax_2014/qtltools/output/array/monocyte_naive.permuted.txt.gz") 
FF_LPS24 = eQTLUtils::importQTLtoolsTable("processed/Fairfax_2014/qtltools/output/array/monocyte_LPS24.permuted.txt.gz") 
FF_IFN = eQTLUtils::importQTLtoolsTable("processed/Fairfax_2014/qtltools/output/array/monocyte_IFN24.permuted.txt.gz") 

#import variant info
gds_file = "results/genotypes/Fairfax_2014/Fairfax_2014_GRCh38.filtered.renamed.gds"
var_info = importVariantInformationFromGDS(gds_file)

#Selected QTL
qtl = dplyr::filter(FF_LPS, group_id == "ENSG00000166750")
fairfax_se = readRDS("results/SummarizedExperiments/array_norm/Fairfax_2014.rds")
assays(fairfax_se) = assays(fairfax_se)[c(1,3)]

#Extract data
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, fairfax_se, "norm_exprs")
selected_genotype = extractVariantGenotypeFromGDS(selected_snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = condition, gene_name = "SLFN5")
qtl_plot = plotQtlRow(selected_data)
ggsave("results/figures/SLFN5_Fairfax_2014.pdf",qtl_plot, width = 4, height = 4)


#PTK2B
qtl = dplyr::filter(FF_naive, group_id == "ENSG00000120899")
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, fairfax_se, "norm_exprs")
selected_genotype = extractVariantGenotypeFromGDS(qtl$snp_id, var_info, gds_file)
selected_data = dplyr::left_join(selected_phenotype, selected_genotype, by = "genotype_id") %>% 
  dplyr::mutate(norm_exp = phenotype_value, genotype_text = factor(genotype_value), 
                condition_name = condition, gene_name = "PTK2B")
qtl_plot = plotQtlRow(selected_data)

#Try to perform finemappign with susieR
qtl = dplyr::filter(FF_naive, group_id == "ENSG00000120899")
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, fairfax_se, "norm_exprs")

#Import genotypes from the region
start = 27337604-100000
end = 27337604+100000
chromosome = "8"
gt = extractGenotypeMatrixFromGDS(chromosome, start, end, var_info, gds_file)

#Map QTLs in the naive condition
naive_data = dplyr::filter(selected_phenotype, condition == "naive")
naive_mat = as.matrix(dplyr::select(naive_data, phenotype_value))
colnames(naive_mat) = "PTK2B"
rownames(naive_mat) = naive_data$genotype_id
gene_pos = dplyr::data_frame(geneid = c("PTK2B"), chr = c("8"), left = c(27337604), right = c(27337604))
snpspos = dplyr::filter(var_info, snp_id %in% rownames(gt)) %>% dplyr::transmute(snpid = snp_id, chr = chromosome, pos)

geno = gt[,rownames(naive_mat)]
naive_res = runMatrixEQTL(t(naive_mat), geno, as.data.frame(snpspos), as.data.frame(gene_pos), cisDist = 1e6, pvOutputThreshold = 1)
qtl_res = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))
ggplot(qtl_res, aes(x = pos, y = -log(pvalue, 10))) + geom_point()

#Use the second lead variant as covariate
covariate = geno["chr8_27371545_T_TTATA",,drop=FALSE]
naive_res = runMatrixEQTL(t(naive_mat), geno, as.data.frame(snpspos), as.data.frame(gene_pos), 
                          covariates = covariate, cisDist = 1e6, pvOutputThreshold = 1)
qtl_res2 = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))
ggplot(qtl_res2, aes(x = pos, y = -log(pvalue, 10))) + geom_point()

dplyr::filter(qtl_res, snps == "chr8_27369273_A_C")
dplyr::filter(qtl_res2, snps == "chr8_27369273_A_C")


covariate = geno["chr8_27369273_A_C",,drop=FALSE]
naive_res = runMatrixEQTL(t(naive_mat), geno, as.data.frame(snpspos), as.data.frame(gene_pos), 
                          covariates = covariate, cisDist = 1e6, pvOutputThreshold = 1)
qtl_res2 = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))
ggplot(qtl_res2, aes(x = pos, y = -log(pvalue, 10))) + geom_point()

dplyr::filter(qtl_res, snps == "chr8_27371545_T_TTATA")
dplyr::filter(qtl_res2, snps == "chr8_27371545_T_TTATA")


#Map QTLs in the LPS24 condition
lps24_data = dplyr::filter(selected_phenotype, condition == "LPS24")
lps24_mat = as.matrix(dplyr::select(lps24_data, phenotype_value))
colnames(lps24_mat) = "PTK2B"
rownames(lps24_mat) = lps24_data$genotype_id
gene_pos = dplyr::data_frame(geneid = c("PTK2B"), chr = c("8"), left = c(27337604), right = c(27337604))
snpspos = dplyr::filter(var_info, snp_id %in% rownames(gt)) %>% dplyr::transmute(snpid = snp_id, chr = chromosome, pos)

lps24_geno = gt[,rownames(lps24_mat)]
naive_res = runMatrixEQTL(t(lps24_mat), lps24_geno, as.data.frame(snpspos), as.data.frame(gene_pos), cisDist = 1e6, pvOutputThreshold = 1)
qtl_res_lps24 = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))
ggplot(qtl_res_lps24, aes(x = pos, y = -log(pvalue, 10))) + geom_point()


covariate = lps24_geno["chr8_27371545_T_TTATA",,drop=FALSE]
naive_res = runMatrixEQTL(t(lps24_mat), lps24_geno, as.data.frame(snpspos), as.data.frame(gene_pos), covariates = covariate, cisDist = 1e6, pvOutputThreshold = 1)
qtl_res_lps24_2 = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))
ggplot(qtl_res_lps24_2, aes(x = pos, y = -log(pvalue, 10))) + geom_point()

dplyr::filter(qtl_res_lps24, snps == "chr8_27369273_A_C")
dplyr::filter(qtl_res_lps24_2, snps == "chr8_27369273_A_C")


covariate = lps24_geno["chr8_27369273_A_C",,drop=FALSE]
naive_res = runMatrixEQTL(t(lps24_mat), lps24_geno, as.data.frame(snpspos), as.data.frame(gene_pos), covariates = covariate, cisDist = 1e6, pvOutputThreshold = 1)
qtl_res_lps24_2 = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))
ggplot(qtl_res_lps24_2, aes(x = pos, y = -log(pvalue, 10))) + geom_point()

dplyr::filter(qtl_res_lps24, snps == "chr8_27371545_T_TTATA")
dplyr::filter(qtl_res_lps24_2, snps == "chr8_27371545_T_TTATA")


#Perform finemapping with susie in both conditions
#Naive
naive_geno_std = geno - apply(geno, 1, mean)
#naive_geno_std = naive_geno_std/apply(naive_geno_std, 1, sd)
naive_geno_std = t(naive_geno_std)

naive_values = as.vector(naive_mat[,1])
naive_values = (naive_values-mean(naive_values))/sd(naive_values)

naive_fitted <- susieR::susie(naive_geno_std, naive_values,
                L = 10,
                estimate_residual_variance = TRUE, 
                estimate_prior_variance = FALSE,
                scaled_prior_variance = 0.1,
                verbose = TRUE)
print(naive_fitted$sets)
naive_cs = naive_fitted$sets$cs$L1
naive_cs_names = colnames(naive_geno_std[,naive_cs])
susieR::susie_plot(naive_fitted, y="PIP")

#LPS24
lps24_geno_std = lps24_geno - apply(lps24_geno, 1, mean)
lps24_geno_std = t(lps24_geno_std)

lps24_values = as.vector(lps24_mat[,1])
lps24_values = (lps24_values-mean(lps24_values))/sd(lps24_values)

lps24_fitted <- susie(lps24_geno_std, lps24_values,
                      L = 10,
                      estimate_residual_variance = TRUE, 
                      estimate_prior_variance = FALSE,
                      scaled_prior_variance = 0.1,
                      verbose = TRUE)
print(lps24_fitted$sets)
lps24_cs = lps24_fitted$sets$cs$L1
lps24_cs_names = colnames(lps24_geno_std[,lps24_cs])
susie_plot(lps24_fitted, y="PIP")

#Combine data and make plots
qtl_df = dplyr::bind_rows(dplyr::mutate(dplyr::arrange(qtl_res, pos), condition = "Control", PIP = naive_fitted$pip), 
                          dplyr::mutate(dplyr::arrange(qtl_res_lps24, pos), condition = "LPS24", PIP = lps24_fitted$pip)) %>%
  dplyr::mutate(credible_set = ifelse(snps %in% naive_cs_names, "CS1",NA)) %>%
  dplyr::mutate(credible_set = ifelse(snps %in% lps24_cs_names, "CS2", credible_set))
plot = ggplot(qtl_df, aes(x = pos, y = -log(pvalue, 10))) + geom_point() + 
  facet_wrap(~condition, ncol = 1) +
  theme_light() +
  xlab("Chromosome 8 position") +
  ylab("-log10 p-value")
ggsave("~/Google Drive/Grant applications/PUT_2019/PTK2B_manhattan.pdf", plot = plot, width = 3.5, height = 3)

#Plot PIP
pip_plot = ggplot(qtl_df, aes(x = pos, y = PIP, colour = credible_set)) + geom_point() + 
  facet_wrap(~condition, ncol = 1) +
  theme_light() +
  xlab("Chromosome 8 position") +
  ylab("Posterior inclusion probability")
ggsave("~/Google Drive/Grant applications/PUT_2019/PTK2B_PIP.pdf", plot = pip_plot, width = 4.5, height = 3)


#Visualise effect sizes
dplyr::mutate(z_score = abs(qnorm(p_nominal/2, mean = 0, sd = 1))*sign(beta)) %>%
  dplyr::mutate(se = beta/abs(z_score))


naive_qtl_effect_sizes = dplyr::bind_rows(dplyr::filter(qtl_res, snps == "chr8_27369273_A_C"),
                 dplyr::filter(qtl_res_lps24, snps == "chr8_27369273_A_C"),
                 dplyr::filter(qtl_res2, snps == "chr8_27369273_A_C"),
                 dplyr::filter(qtl_res_lps24_2, snps == "chr8_27369273_A_C")) %>%
  dplyr::mutate(analysis_type = c("marginal","marginal","adjusted", "adjusted")) %>%
  dplyr::mutate(condition = c("Control","LPS24","Control", "LPS24")) %>%
  dplyr::mutate(z_score = abs(qnorm(pvalue/2, mean = 0, sd = 1))*sign(beta)) %>%
  dplyr::mutate(se = beta/abs(z_score))

lps24_qtl_effect_sizes = dplyr::bind_rows(dplyr::filter(qtl_res, snps == "chr8_27371545_T_TTATA"),
                 dplyr::filter(qtl_res_lps24, snps == "chr8_27371545_T_TTATA"),
                 dplyr::filter(qtl_res2, snps == "chr8_27371545_T_TTATA"),
                 dplyr::filter(qtl_res_lps24_2, snps == "chr8_27371545_T_TTATA")) %>%
  dplyr::mutate(analysis_type = c("marginal","marginal","adjusted", "adjusted")) %>%
  dplyr::mutate(condition = c("Control","LPS24","Control", "LPS24")) %>%
  dplyr::mutate(z_score = abs(qnorm(pvalue/2, mean = 0, sd = 1))*sign(beta)) %>%
  dplyr::mutate(se = beta/abs(z_score))

table = dplyr::bind_rows(naive_qtl_effect_sizes, lps24_qtl_effect_sizes) %>% 
  dplyr::mutate(analysis_type = ifelse(analysis_type == "adjusted", "conditional", analysis_type)) %>%
  dplyr::mutate(analysis_type = factor(analysis_type, levels = c("marginal", "conditional"))) %>%
  dplyr::mutate(lead_var = ifelse(snps == "chr8_27369273_A_C", "Control eQTL", "LPS24 eQTL"))
effect_plot = ggplot(table, aes(x = condition, y = beta, ymax = beta+se, ymin = beta-se)) + 
  geom_point() + 
  geom_errorbar(width = 0.1) + 
  geom_hline(yintercept = 0) + 
  facet_grid(lead_var~analysis_type) + 
  ylab("Effect on PTK2B expression") +
  coord_flip()
ggsave("~/Google Drive/Grant applications/PUT_2019/PTK2B_effect_estimates.pdf", effect_plot, width = 3, height = 2.5)

#Do the same analysis for CD40
qtl = dplyr::filter(FF_naive, group_id == "ENSG00000101017")
qtl$p_nominal

#Import genotypes from the region
start = qtl$snp_start-100000
end = qtl$snp_start+100000
chromosome = qtl$snp_chr
gt = extractGenotypeMatrixFromGDS(chromosome, start, end, var_info, gds_file)
selected_phenotype = extractPhentypeFromSE(qtl$phenotype_id, fairfax_se, "norm_exprs")

#Map QTLs in the naive condition
naive_data = dplyr::filter(selected_phenotype, condition == "naive")
naive_mat = as.matrix(dplyr::select(naive_data, phenotype_value))
colnames(naive_mat) = "CD40"
rownames(naive_mat) = naive_data$genotype_id
gene_pos = dplyr::data_frame(geneid = c("CD40"), chr = chromosome, left = qtl$snp_start, right = qtl$snp_start)
snpspos = dplyr::filter(var_info, snp_id %in% rownames(gt)) %>% dplyr::transmute(snpid = snp_id, chr = chromosome, pos)

#Extract genotype data
geno = gt[,rownames(naive_mat)]
geno_std = geno - apply(geno, 1, mean)
geno_std = t(geno_std)

exprs = naive_mat[,1]
std_exprs = (exprs-mean(exprs))/sd(exprs)

fitted <- susie(geno_std, std_exprs,
                      L = 10,
                      estimate_residual_variance = TRUE, 
                      estimate_prior_variance = FALSE,
                      scaled_prior_variance = 0.1,
                      verbose = TRUE)
print(fitted$sets)
susie_plot(fitted, y="PIP")

naive_res = runMatrixEQTL(t(naive_mat), geno, as.data.frame(snpspos), as.data.frame(gene_pos), cisDist = 1e6, pvOutputThreshold = 1)
qtl_res = naive_res$cis$eqtls %>% dplyr::as_tibble() %>%
  dplyr::left_join(snpspos, by = c("snps" = "snpid"))

cs1 = rownames(geno)[fitted$sets$cs$L1]
cs2 = rownames(geno)[fitted$sets$cs$L2]
qtl_res = dplyr::mutate(dplyr::arrange(qtl_res, pos), PIP = fitted$pip) %>%
  dplyr::mutate(credible_set = ifelse(snps %in% cs1, "CS1", NA)) %>%
  dplyr::mutate(credible_set = ifelse(snps %in% cs2, "CS2", credible_set))

association_plot = ggplot(qtl_res, aes(x = pos, y = -log(pvalue, 10))) + 
  geom_point() + 
  theme_light() + 
  xlab("Chromosome 20 position") +
  ylab("-log10 p-value")
ggsave("~/Google Drive/Grant applications/PUT_2019/CD40_pvalues.pdf", association_plot, width = 4, height = 2.5)

pip_plot = ggplot(qtl_res, aes(x = pos, y = PIP, color = credible_set)) + 
  geom_point() + 
  theme_light() + 
  xlab("Chromosome 20 position") + 
  ylab("Posterior inclusion probability")
ggsave("~/Google Drive/Grant applications/PUT_2019/CD40_PIP.pdf", pip_plot, width = 5, height = 2.5)




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

