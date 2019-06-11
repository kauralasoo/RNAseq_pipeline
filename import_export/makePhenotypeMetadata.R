library("data.table")
library("devtools")
load_all("../eQTLUtils/")


#### Gene counts (featureCounts) ####
#Specify required phenotype metadata columns
required_phenotype_meta_columns = c("phenotype_id","quant_id","group_id","gene_id","chromosome","gene_start",
                               "gene_end","strand","gene_name","gene_type","gene_version","phenotype_pos")
required_gene_meta_columns = c(required_phenotype_meta_columns, "phenotype_gc_content", "phenotype_length")

#Import gene metadata
transcript_meta = readr::read_tsv("annotations/eQTLCatalogue/v0.1/Homo_sapiens.GRCh38.96_biomart_download.txt.gz", col_types = "ccccccciciiciiiiccccccccidccccii")
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


