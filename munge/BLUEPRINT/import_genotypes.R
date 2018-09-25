library("SNPRelate")
library("GDSArray")

#Make a GDS file
snpgdsVCF2GDS("results/genotypes/BLUEPRINT/BLUEPRINT_06092016_GRCh38.sorted.MAF05.vcf.gz", 
              "results/genotypes/BLUEPRINT/BLUEPRINT_06092016_GRCh38.sorted.MAF05.gds", 
              method = "copy.num.of.ref")

a = GDSFile("results/genotypes/BLUEPRINT/BLUEPRINT_06092016_GRCh38.sorted.MAF05.gds")

geno = GDSArray("results/genotypes/BLUEPRINT/BLUEPRINT_06092016_GRCh38.sorted.MAF05.gds", "genotype")


GDSArray("results/genotypes/BLUEPRINT/BLUEPRINT_06092016_GRCh38.sorted.MAF05.gds", "snp.chromosome")
