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
ggplot(pca_df, aes(x = PC1, y = PC2)) + geom_point()

#GENCORD
pca_df = importQTLtoolsPCA("processed/GENCORD/qtltools/input/vcf/LCL.geno.pca")
ggplot(pca_df, aes(x = PC1, y = PC2)) + geom_point()
