library("readr")
library("dplyr")
library("devtools")
library("data.table")
load_all("../eQTLUtils/")

#Import biomart data
transcript_meta = importBiomartMetadata("annotations/eQTLCatalogue/v0.1/Homo_sapiens.GRCh38.96_biomart_download.txt.gz")

#Import Leafcutter count matrix
lc_matrix = utils::read.csv("results/Schwartzentruber_2018/leafcutter_perind_numers.counts.formatted.gz", sep = "", stringsAsFactors = FALSE)
leafcutter_meta = leafcutterAnnotateIntrons(lc_matrix$phenotype_id, 
        intron_annotation_path = "annotations/eQTLCatalogue/v0.1/leafcutter_annotations/gencode_v30_all_introns.bed.gz",
        transcript_meta)

#Export Leafcutter gene metadata file
write.table(leafcutter_meta, "results/Schwartzentruber_2018/leafcutter_cluster_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Import sample metadata
sample_meta = read.table("../SampleArcheology/studies/cleaned/Schwartzentruber_2018.tsv", stringsAsFactors = F, header = T, sep = "\t") %>% 
  dplyr::as_tibble() %>% dplyr::filter(rna_qc_passed, genotype_qc_passed)

#Make Summarized experiment
lc_se = makeSummarizedExperimentFromCountMatrix(lc_matrix, leafcutter_meta, sample_meta, assay_name = "counts", quant_method = "leafcutter")

#Normalise
lc_norm = qtltoolsPrepareSE(lc_se[,1:10], quant_method = "leafcutter")

#Write count matrix as tab-separated file
gz1 = gzfile("results/Schwartzentruber_2018/lc_norm.tsv.gz","w") 
write.table(assays(lc_norm)$qnorm, gz1, sep = "\t", quote = F)
close(gz1)

