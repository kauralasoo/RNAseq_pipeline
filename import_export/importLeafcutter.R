library("readr")
library("dplyr")
library("devtools")
library("data.table")
load_all("../eQTLUtils/")

#Import biomart data
transcript_meta = importBiomartMetadata("annotations/Ensembl94_biomart_download.txt.gz")

#Import Leafcutter count matrix
geuvadis_lc = read.table("processed/GEUVADIS/leafcutter/leafcutter_perind_numers.counts.gz", sep = " ")
leafcutter_meta = leafcutterAnnotateIntrons(rownames(geuvadis_lc), 
        intron_annotation_path = "annotations/gencode_v29_leafcutter_annot_all_introns.no_chr.bed.gz",
        transcript_meta)

#Import sample metadata
geuvadis_meta = read.table("../SampleArcheology/studies/cleaned/GEUVADIS_EUR.tsv", stringsAsFactors = F, header = T) %>% 
  dplyr::as_tibble()

#Make Summarized experiment
geuvadis_lc_se = makeSummarizedExperiment(geuvadis_lc, leafcutter_meta, geuvadis_meta, assay_name = "counts")

#Normalise
geuvadis_norm = qtltoolsPrepareSE(geuvadis_lc_se[,1:10], quant_method = "LeafCutter")

#Export Leafcutter gene metadata file
write.table(leafcutter_meta, "metadata/gene_metadata/GEUVADIS_leafcutter_cluster_metadata.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#Write count matrix as tab-separated file
gz1 = gzfile("results/expression_matrices/LeafCutter/GEUVADIS.tsv.gz","w") 
write.table(geuvadis_lc, gz1, sep = "\t", quote = F)
close(gz1)

  
