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
geuvadis_meta = read.table("metadata/cleaned/GEUVADIS.tsv", stringsAsFactors = F, header = T) %>% 
  dplyr::as_tibble()

geuvadis_lc_se = makeSummarizedExperiment(geuvadis_lc, leafcutter_meta, geuvadis_meta, assay_name = "counts")




  



