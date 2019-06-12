library("data.table")
library("devtools")
load_all("../eQTLUtils/")

#### Gene counts (featureCounts) ####
#Specify required phenotype metadata columns
required_phenotype_meta_columns = c("phenotype_id","quant_id","group_id","gene_id","chromosome","gene_start",
                               "gene_end","strand","gene_name","gene_type","gene_version","phenotype_pos")
required_gene_meta_columns = c(required_phenotype_meta_columns, "phenotype_gc_content", "phenotype_length")

#Import gene metadata
transcript_meta = eQTLUtils::importBiomartMetadata("annotations/eQTLCatalogue/v0.1/Homo_sapiens.GRCh38.96_biomart_download.txt.gz")
gene_metadata = extractGeneMetadataFromBiomartFile(transcript_meta)

#Import featureCounts gene lengths
lengths = read.table("annotations/helper_files/SRR4295627.sorted_gene.featureCounts.txt", 
                     header = T, stringsAsFactors = F) %>% 
  dplyr::as_tibble() %>%
  dplyr::transmute(phenotype_id = Geneid, phenotype_length = Length) %>%
  dplyr::filter(!(phenotype_id %like% "PAR_Y")) %>%
  reformatPhenotypeId()

#Make gene metadata
final_gene_metadata = dplyr::left_join(lengths, gene_metadata, by = "phenotype_id") %>%
  dplyr::rename(phenotype_gc_content = gene_gc_content) %>%
  dplyr::select(required_gene_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("metadata/phenotype_metadata/gene_counts_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(final_gene_metadata, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

#### Exon counts ####
exon_nuc_content = readr::read_tsv("annotations/eQTLCatalogue/v0.1/DEXSeq/gencode.v30.annotation.no_chr.patched_contigs.DEXSeq.nuc_content.txt",
                                   col_types = "ccciiccccddiiiiiii")
nuc_content = exon_nuc_content[,c(1,4,5,11)]
colnames(nuc_content) = c("seqnames", "start", "end", "phenotype_gc_content")
nuc_content = dplyr::distinct(nuc_content)

#Subset of gene metadata
gene_metadata_subset = dplyr::select(gene_metadata, gene_id, chromosome, gene_start, gene_end, 
                                     strand, gene_name, gene_type, gene_version)
  
#Make exon metadata
exon_metadata = rtracklayer::import.gff("annotations/eQTLCatalogue/v0.1/DEXSeq/gencode.v30.annotation.no_chr.patched_contigs.DEXSeq.gff") %>%
  as.data.frame() %>% dplyr::as_tibble() %>%
  dplyr::filter(type != "aggregate_gene") %>%
  dplyr::mutate(phenotype_id = paste(gene_id, seqnames, start, end, sep = "_")) %>%
  dplyr::mutate(phenotype_pos = floor((start + end)/2)) %>%
  dplyr::transmute(phenotype_id, gene_id, seqnames, start, end, phenotype_length = width, phenotype_pos) %>%
  dplyr::mutate(seqnames = as.character(seqnames)) %>%
  dplyr::left_join(nuc_content, by = c("seqnames", "start", "end")) %>%
  dplyr::filter(!grepl("+", gene_id, fixed = TRUE)) %>%
  removeGeneVersion() %>%
  dplyr::mutate(quant_id = gene_id, group_id = gene_id) %>%
  dplyr::left_join(gene_metadata_subset, by = "gene_id") %>%
  dplyr::select(-seqnames, -start, -end) %>%
  dplyr::select(required_gene_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("metadata/phenotype_metadata/exon_counts_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(exon_metadata, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)


#### GENCODE transcripts ####
tx_gene_map = dplyr::transmute(transcript_meta, phenotype_id = transcript_id, gene_id)
tx_estimates = readr::read_tsv("annotations/helper_files/gencode.v30.transcripts_TPM_merged.txt") %>% reformatPhenotypeId()
gencode_transcripts_meta = dplyr::select(tx_estimates, phenotype_id) %>% 
  dplyr::left_join(tx_gene_map) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("metadata/phenotype_metadata/transcript_usage_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(gencode_transcripts_meta, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

### txrevise events ####
list = as.list(c("annotations/helper_files/merged_counts/txrevise.grp_1.contained/txrevise.grp_1.contained_TPM_merged.txt", 
                 "annotations/helper_files/merged_counts/txrevise.grp_2.contained/txrevise.grp_2.contained_TPM_merged.txt",
                 "annotations/helper_files/merged_counts/txrevise.grp_1.upstream/txrevise.grp_1.upstream_TPM_merged.txt",
                 "annotations/helper_files/merged_counts/txrevise.grp_2.upstream/txrevise.grp_2.upstream_TPM_merged.txt",
                 "annotations/helper_files/merged_counts/txrevise.grp_1.downstream/txrevise.grp_1.downstream_TPM_merged.txt",
                 "annotations/helper_files/merged_counts/txrevise.grp_2.downstream/txrevise.grp_2.downstream_TPM_merged.txt"))
event_quants = purrr::map_df(list, ~readr::read_tsv(.))

#Extract relevant infromation from the txrevise events
txrevise_meta = dplyr::select(event_quants, phenotype_id) %>% 
  tidyr::separate(phenotype_id, c("gene_id", "grp", "pos", "transcript"), sep = "\\.", remove = F) %>%
  dplyr::mutate(quant_id = paste(gene_id, grp, pos, sep = "."), group_id = paste(gene_id, pos, sep = ".")) %>% 
  dplyr::select(phenotype_id, quant_id, group_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content, -group_id, -quant_id), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything())

#Save expression matrix
gz2 = gzfile("metadata/phenotype_metadata/txrevise_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(txrevise_meta, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

#### Human HT-12 V4 metadata ####
human_ht12 = readr::read_tsv("metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz") %>%
  dplyr::select(phenotype_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything()) %>%
  dplyr::filter(!is.na(gene_name))

#Save metadata file
gz2 = gzfile("metadata/phenotype_metadata/HumanHT-12_V4_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(human_ht12, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)

#### Human HT-12 V4 metadata ####
affy = readr::read_tsv("metadata/gene_metadata/Affy_Human_Gene_1_0_ST_gene_metadata.txt.gz", col_types = "ccccciiiccddi") %>%
  dplyr::select(phenotype_id, gene_id) %>%
  dplyr::left_join(dplyr::select(gene_metadata, -phenotype_id, -gene_gc_content), by = "gene_id") %>%
  dplyr::select(required_phenotype_meta_columns, dplyr::everything()) %>%
  dplyr::filter(!is.na(gene_name))

#Save metadata file
gz2 = gzfile("metadata/phenotype_metadata/Affy_Human_Gene_1_0_ST_Ensembl_96_phenotype_metadata.tsv.gz", "w")
write.table(affy, gz2, sep = "\t", quote = FALSE, row.names = F)
close(gz2)



