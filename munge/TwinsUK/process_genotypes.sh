# Add INFO scores and standard variant IDs
bcftools +impute-info TwinsUK_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o TwinsUK_GRCh38.annotated.vcf.gz

bcftools filter -i 'MAF[0] > 0.05' TwinsUK_GRCh38.annotated.vcf.gz -Oz -o TwinsUK_GRCh38.filtered.vcf.gz