#Merge all chormosomes
bcftools concat GEUVADIS_445_samples.chr1.vcf.gz GEUVADIS_445_samples.chr2.vcf.gz GEUVADIS_445_samples.chr10.vcf.gz GEUVADIS_445_samples.chr11.vcf.gz GEUVADIS_445_samples.chr12.vcf.gz GEUVADIS_445_samples.chr13.vcf.gz GEUVADIS_445_samples.chr14.vcf.gz GEUVADIS_445_samples.chr15.vcf.gz GEUVADIS_445_samples.chr16.vcf.gz GEUVADIS_445_samples.chr17.vcf.gz GEUVADIS_445_samples.chr18.vcf.gz GEUVADIS_445_samples.chr19.vcf.gz GEUVADIS_445_samples.chr20.vcf.gz GEUVADIS_445_samples.chr21.vcf.gz GEUVADIS_445_samples.chr22.vcf.gz GEUVADIS_445_samples.chr3.vcf.gz GEUVADIS_445_samples.chr4.vcf.gz GEUVADIS_445_samples.chr5.vcf.gz GEUVADIS_445_samples.chr6.vcf.gz GEUVADIS_445_samples.chr7.vcf.gz GEUVADIS_445_samples.chr8.vcf.gz GEUVADIS_445_samples.chr9.vcf.gz GEUVADIS_445_samples.chrX.vcf.gz -O z -o GEUVADIS_445_samples.merged.vcf.gz

#Split multi-allelic variants
bcftools norm -m-any GEUVADIS_445_samples.merged.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o GEUVADIS_GRCh38_filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags GEUVADIS_GRCh38_filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\tNA\n' | gzip > GEUVADIS_GRCh38.variant_information.txt.gz
