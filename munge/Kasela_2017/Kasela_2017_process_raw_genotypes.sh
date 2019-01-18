#Merge X in the PAR region
plink --file CTG_omniSNP_unik313ind --merge-x --make-bed --out Kasela_2017_GRCh37

#Update build
update_build.sh Kasela_2017_GRCh37 HumanOmniExpress-12v1_A-b37.strand Kasela_2017_GRCh37_strand

#Convert to VCF
plink --bfile Kasela_2017_GRCh37_strand --recode vcf-iid --out Kasela_2017_GRCh37

#Rename 23 and 24 to X and Y
bcftools annotate --rename-chrs ~/hpc/projects/RNAseq_pipeline/genotype_scripts/chromsome_names.txt Kasela_2017_GRCh37.vcf -Oz -o Kasela_2017_GRCh37_renamed.vcf.gz

#Fix the ref allele
bcftools index Kasela_2017_GRCh37_renamed.vcf.gz
bcftools +fixref Kasela_2017_GRCh37_renamed.vcf.gz -Oz -o Kasela_2017_GRCh37_renamed_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort the vcf file 
bcftools sort Kasela_2017_GRCh37_renamed_fixref.vcf.gz -Oz -o Kasela_2017_GRCh37_renamed_fixref_sorted.vcf.gz

#Remove all variants with missing ALT alleles
bcftools filter -e 'ALT="."' Kasela_2017_GRCh37_renamed_fixref_sorted.vcf.gz -Oz -o Kasela_2017_GRCh37_renamed_fixref_sorted_noALT.vcf.gz

#Add all INFO tags
bcftools +fill-tags Kasela_2017_GRCh37_renamed_fixref_sorted_noALT.vcf.gz -Oz -o Kasela_2017_GRCh37_renamed_fixref_sorted_noALT_tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Kasela_2017_GRCh37_renamed_fixref_sorted_noALT_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o Kasela_2017_GRCh37_filtered.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa Kasela_2017_GRCh37_filtered.vcf.gz -Oz -o Kasela_2017_GRCh37_filtered_noref.vcf.gz

#Remove duplicates and multi-allelics
bcftools norm -d all Kasela_2017_GRCh37_filtered_noref.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o Kasela_2017_GRCh37_genotyped.vcf.gz

#Extract chromosomes
bcftools index Kasela_2017_GRCh37_genotyped.vcf.gz
bcftools view -r 1 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr1.vcf.gz
bcftools view -r 2 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr2.vcf.gz
bcftools view -r 3 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr3.vcf.gz
bcftools view -r 4 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr4.vcf.gz
bcftools view -r 5 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr5.vcf.gz
bcftools view -r 6 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr6.vcf.gz
bcftools view -r 7 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr7.vcf.gz
bcftools view -r 8 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr8.vcf.gz
bcftools view -r 9 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr9.vcf.gz
bcftools view -r 10 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr10.vcf.gz
bcftools view -r 11 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr11.vcf.gz
bcftools view -r 12 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr12.vcf.gz
bcftools view -r 13 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr13.vcf.gz
bcftools view -r 14 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr14.vcf.gz
bcftools view -r 15 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr15.vcf.gz
bcftools view -r 16 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr16.vcf.gz
bcftools view -r 17 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr17.vcf.gz
bcftools view -r 18 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr18.vcf.gz
bcftools view -r 19 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr19.vcf.gz
bcftools view -r 20 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr20.vcf.gz
bcftools view -r 21 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr21.vcf.gz
bcftools view -r 22 Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chr22.vcf.gz
bcftools view -r X Kasela_2017_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Kasela_2017_GRCh37_chrX.vcf.gz