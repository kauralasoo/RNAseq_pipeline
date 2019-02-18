# Add INFO scores and standard variant IDs
bcftools +impute-info TwinsUK_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o TwinsUK_GRCh38.annotated.vcf.gz

#Keep only unrelated individuals
bcftools view -S ~/rocket/projects/GEUVADIS/RNAseq_pipeline/metadata/TwinsUK/TwinsUK_unrelated_individuals.txt TwinsUK_GRCh38.annotated.vcf.gz | bcftools filter -i 'MAF[0] > 0.01' -Oz -o TwinsUK_GRCh38.filtered.vcf.gz
bcftools index TwinsUK_GRCh38.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags TwinsUK_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\tNA\n' | gzip > TwinsUK_GRCh38.variant_information.txt.gz
