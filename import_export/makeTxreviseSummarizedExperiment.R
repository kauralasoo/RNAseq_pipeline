library("tximport")
library("dplyr")
library("devtools")
library("tidyr")
load_all("../seqUtils/")

#### txrevise ####
sample_metadata = read.table("../RNAseq_pipeline/metadata/cleaned/GEUVADIS.tsv", header = TRUE, stringsAsFactors = F)

#Import txrevise quant results
annotations = c("txrevise.grp_1.upstream","txrevise.grp_2.upstream",
                "txrevise.grp_1.contained","txrevise.grp_2.contained",
                "txrevise.grp_1.downstream","txrevise.grp_2.downstream")
file_names = paste0("../RNAseq_pipeline/processed/GEUVADIS/matrices/", annotations, ".salmon_txrevise.rds")
file_list = setNames(as.list(file_names), annotations)
quant_matrix = purrr::map(file_list, readRDS) %>%
  purrr::map(~.$abundance) %>%
  purrr::reduce(rbind)

#Construct event metadata
transcript_meta = importBiomartMetadata("../RNAseq_pipeline/annotations/Ensembl92_biomart_download.txt.gz")
row_data = constructTxreviseRowData(rownames(quant_matrix), transcript_meta)

#Make SummarizedExperiment
txrevise_se = makeSummarizedExperiment(quant_matrix, row_data, sample_metadata, assay_name = "tpms")

#Calculate transcript ratios
txrevise_usage = normaliseSE_ratios(txrevise_se, assay_name = "tpms")
saveRDS(txrevise_usage, "results/SummarizedExperiments/GEUVADIS_txrevise.rds")








