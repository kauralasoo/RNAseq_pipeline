#Extract macrophage samples from the large HipSci vcf file
bcftools view -S donor_ids.txt ~/hpc/datasets/controlled_access/HipSci/genotypes/Michigan_GRCh37_1KGPhase3_081018/GRCh38/HipSci_GRCh38.filtered.vcf.gz -Oz -o Alasoo_2018.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags Alasoo_2018_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\t%R2\n' | gzip > Alasoo_2018_GRCh38.variant_information.txt.gz
