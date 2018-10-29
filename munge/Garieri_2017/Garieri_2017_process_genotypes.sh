#Merge all chormosomes
bcftools concat Garieri_1000G.chr1.vcf.gz Garieri_1000G.chr2.vcf.gz Garieri_1000G.chr3.vcf.gz Garieri_1000G.chr4.vcf.gz Garieri_1000G.chr5.vcf.gz Garieri_1000G.chr6.vcf.gz Garieri_1000G.chr7.vcf.gz Garieri_1000G.chr8.vcf.gz Garieri_1000G.chr9.vcf.gz Garieri_1000G.chr10.vcf.gz Garieri_1000G.chr11.vcf.gz Garieri_1000G.chr12.vcf.gz Garieri_1000G.chr13.vcf.gz Garieri_1000G.chr14.vcf.gz Garieri_1000G.chr15.vcf.gz Garieri_1000G.chr16.vcf.gz Garieri_1000G.chr17.vcf.gz Garieri_1000G.chr18.vcf.gz Garieri_1000G.chr19.vcf.gz Garieri_1000G.chr20.vcf.gz Garieri_1000G.chr21.vcf.gz Garieri_1000G.chr22.vcf.gz -Oz -o Garieri_1000G.vcf.gz
bcftools index Garieri_1000G.vcf.gz

#Split multi-allelic variants
bcftools norm -m-any Garieri_1000G.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Garieri_1000G_renamed.vcf.gz

#Remove all INFO tags
bcftools annotate -x INFO Garieri_1000G_renamed.vcf.gz -Oz -o Garieri_1000G_renamed_no_tags.vcf.gz

#Merge the two files
bcftools merge 1000G/Garieri_1000G_renamed_no_tags.vcf.gz GENCORD_subset.vcf.gz -Oz -o Garieri_merged.vcf.gz

#Keep genotypes that are not missing
bcftools filter -i 'F_MISSING < 0.05 & MAF[0] > 0.01' Garieri_merged.vcf.gz | bcftools annotate -x FORMAT/GP | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Garieri_filtered.vcf.gz