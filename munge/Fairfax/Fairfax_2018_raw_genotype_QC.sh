#Rename chromosmes to GRCh37
bcftools annotate --rename-chrs chr_names2.txt out.vcf -Oz -o Fairfax_GRCh37.vcf.gz

#Sort
bcftools sort Fairfax_GRCh37.vcf.gz -Oz -o Fairfax_GRCh37_sorted.vcf.gz
bcftools index Fairfax_GRCh37_sorted.vcf.gz

#Keep chromosomes
bcftools view -r 1,10,11,12,13,14,15,16,17,18,19,2,20,21,22,3,4,5,6,7,8,9,X,Y Fairfax_GRCh37_sorted.vcf.gz -Oz -o Fairfax_GRCh37_filtered.vcf.gz

#Identify individuals with high level of missing genotypes
module load vcftools-0.1.13
vcftools --gzvcf Fairfax_GRCh37_filtered.vcf.gz --missing-indv --out Fairfax_missing

#Add all INFO tags
bcftools +fill-tags Fairfax_GRCh37_filtered.vcf.gz -Oz -o Fairfax_GRCh37_tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Fairfax_GRCh37_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o Fairfax_GRCh37_filtered.vcf.gz

#Fix ref and alt alleles
bcftools index Fairfax_GRCh37_filtered.vcf.gz
bcftools +fixref Fairfax_GRCh37_filtered.vcf.gz -Oz -o Fairfax_GRCh37_filtered_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa Fairfax_GRCh37_filtered_fixref.vcf.gz -Oz -o Fairfax_GRCh37_filtered_nonref.vcf.gz

#Sort
bcftools sort Fairfax_GRCh37_filtered_nonref.vcf.gz -Oz -o Fairfax_GRCh37_genotyped.vcf.gz
bcftools index Fairfax_GRCh37_genotyped.vcf.gz

#Extract chromosomes
bcftools view -r 1 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr1.vcf.gz
bcftools view -r 2 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr2.vcf.gz
bcftools view -r 3 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr3.vcf.gz
bcftools view -r 4 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr4.vcf.gz
bcftools view -r 5 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr5.vcf.gz
bcftools view -r 6 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr6.vcf.gz
bcftools view -r 7 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr7.vcf.gz
bcftools view -r 8 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr8.vcf.gz
bcftools view -r 9 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr9.vcf.gz
bcftools view -r 10 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr10.vcf.gz
bcftools view -r 11 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr11.vcf.gz
bcftools view -r 12 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr12.vcf.gz
bcftools view -r 13 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr13.vcf.gz
bcftools view -r 14 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr14.vcf.gz
bcftools view -r 15 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr15.vcf.gz
bcftools view -r 16 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr16.vcf.gz
bcftools view -r 17 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr17.vcf.gz
bcftools view -r 18 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr18.vcf.gz
bcftools view -r 19 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr19.vcf.gz
bcftools view -r 20 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr20.vcf.gz
bcftools view -r 21 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr21.vcf.gz
bcftools view -r 22 Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chr22.vcf.gz
bcftools view -r X Fairfax_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Fairfax_2018_GRCh37_chrX.vcf.gz