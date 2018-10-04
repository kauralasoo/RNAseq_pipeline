library("readr")
library("dplyr")
library("tidyr")
library("ggplot2")
library("cqn")
library("SummarizedExperiment")
library("ggplot2")
library("data.table")
library("devtools")
load_all("../eQTLUtils/")

#Quach_2016
pca_df = importQTLtoolsPCA("processed/Quach_2016/qtltools/input/featureCounts/vcf/LPS.geno.pca")
ggplot(pca_df, aes(x = PC5, y = PC6)) + geom_point()

#GENCORD
pca_df = importQTLtoolsPCA("processed/GENCORD/qtltools/input/vcf/LCL.geno.pca")
ggplot(pca_df, aes(x = PC1, y = PC2, label = sample_id)) + geom_point() + geom_text()

n = importQTLtoolsTable("processed/Quach_2016/qtltools/output/naive.permuted.txt.gz")
lcl = importQTLtoolsTable("processed/GENCORD/qtltools/output/featureCounts/LCL.permuted.txt.gz")

lcl = importQTLtoolsTable("processed/GEUVADIS/")
