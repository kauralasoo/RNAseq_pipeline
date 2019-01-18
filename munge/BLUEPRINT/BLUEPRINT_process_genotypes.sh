
# Add standard variant IDs
bcftools annotate BLUEPRINT_06092016_GRCh38.sorted.ref.filtered.vcf.gz --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o BLUEPRINT_06092016_GRCh38.annotated.vcf.gz

#Keep only unrelated individuals
bcftools filter BLUEPRINT_06092016_GRCh38.annotated.vcf.gz -i 'MAF[0] > 0.01' -Oz -o BLUEPRINT_06092016_GRCh38_filtered.vcf.gz
bcftools index BLUEPRINT_06092016_GRCh38_filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags BLUEPRINT_06092016_GRCh38_filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\tNA\n' | gzip > BLUEPRINT_06092016_GRCh38.variant_information.txt.gz


