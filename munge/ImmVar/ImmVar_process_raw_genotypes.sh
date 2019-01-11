#Update build
# module load plink-1.9.0
# update_build.sh immvar_omxex_zcall_HG19_Oct2018 HumanOmniExpressExome-8v1_A-b37.Ilmn.strand ImmVar_b37

#Convert to VCF
plink --bfile immvar_omxex_zcall_HG19_Oct2018 --recode vcf-iid --out ImmVar_GRCh37

#Rename 23 and 24 to X and Y
bcftools annotate --rename-chrs ~/hpc/projects/RNAseq_pipeline/genotype_scripts/chromsome_names.txt ImmVar_GRCh37.vcf -Oz -o ImmVar_GRCh37_renamed.vcf.gz

#Count missing variants
module load vcftools-0.1.13
vcftools --vcf ImmVar_GRCh37.vcf --missing-indv --out ImmVar_missing_individuals
vcftools --vcf ImmVar_GRCh37.vcf --missing-site --out ImmVar_missing_sites

#Remove all variants with missing ALT alleles
bcftools filter -e 'ALT="."' ImmVar_GRCh37_renamed.vcf.gz -Oz -o ImmVar_GRCh37_noALT.vcf.gz

#Add all INFO tags
bcftools +fill-tags ImmVar_GRCh37_noALT.vcf.gz -Oz -o ImmVar_GRCh37_tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' ImmVar_GRCh37_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o ImmVar_GRCh37_tags_filtered.vcf.gz

#Fix the ref allele
bcftools index ImmVar_GRCh37_tags_filtered.vcf.gz
bcftools +fixref ImmVar_GRCh37_tags_filtered.vcf.gz -Oz -o ImmVar_GRCh37_tags_filtered_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort the vcf file 
bcftools sort ImmVar_GRCh37_tags_filtered_fixref.vcf.gz -Oz -o ImmVar_GRCh37_tags_filtered_fixref_sorted.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa ImmVar_GRCh37_tags_filtered_fixref_sorted.vcf.gz -Oz -o  ImmVar_GRCh37_tags_filtered_fixref_sorted_noref.vcf.gz

#Remove duplicates and multi-allelics
bcftools norm -d all ImmVar_GRCh37_tags_filtered_fixref_sorted_noref.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o ImmVar_GRCh37_tags_filtered_fixref_sorted_noref_nodup.vcf.gz

#Rename
mv ImmVar_GRCh37_tags_filtered_fixref_sorted_noref_nodup.vcf.gz ImmVar_GRCh37_genotyped.vcf.gz

#Extract chromosomes
mkdir by_chr
bcftools index ImmVar_GRCh37_genotyped.vcf.gz
bcftools view -r 1 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr1.vcf.gz
bcftools view -r 2 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr2.vcf.gz
bcftools view -r 3 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr3.vcf.gz
bcftools view -r 4 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr4.vcf.gz
bcftools view -r 5 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr5.vcf.gz
bcftools view -r 6 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr6.vcf.gz
bcftools view -r 7 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr7.vcf.gz
bcftools view -r 8 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr8.vcf.gz
bcftools view -r 9 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr9.vcf.gz
bcftools view -r 10 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr10.vcf.gz
bcftools view -r 11 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr11.vcf.gz
bcftools view -r 12 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr12.vcf.gz
bcftools view -r 13 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr13.vcf.gz
bcftools view -r 14 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr14.vcf.gz
bcftools view -r 15 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr15.vcf.gz
bcftools view -r 16 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr16.vcf.gz
bcftools view -r 17 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr17.vcf.gz
bcftools view -r 18 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr18.vcf.gz
bcftools view -r 19 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr19.vcf.gz
bcftools view -r 20 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr20.vcf.gz
bcftools view -r 21 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr21.vcf.gz
bcftools view -r 22 ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chr22.vcf.gz
bcftools view -r X ImmVar_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ImmVar_GRCh37_chrX.vcf.gz
