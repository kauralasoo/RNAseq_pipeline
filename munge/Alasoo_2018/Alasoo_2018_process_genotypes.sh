#Extract macrophage samples from the large HipSci vcf file
bcftools view -S donor_ids.txt ~/hpc/datasets/controlled_access/HipSci/genotypes/Michigan_GRCh37_1KGPhase3_081018/GRCh38/HipSci_GRCh38.filtered.vcf.gz -Oz -o Alasoo_2018.filtered.vcf.gz

