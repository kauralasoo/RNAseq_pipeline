# Add INFO scores and standard variant IDs
bcftools +impute-info TwinsUK_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o TwinsUK_GRCh38.annotated.vcf.gz

#Keep only unrelated individuals
bcftools view -S  ~/rocket/projects/GEUVADIS/RNAseq_pipeline/metadata/TwinsUK/TwinsUK_unrelated_individuals.txt TwinsUK_GRCh38.annotated.vcf.gz | bcftools filter -i 'MAF[0] > 0.05' -Oz -o TwinsUK_GRCh38.filtered.vcf.gz