#Add all INFO tags
bcftools +fill-tags HipSci_GRCh37_OA.vcf.gz -Oz -o HipSci_GRCh37_OA.tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' HipSci_GRCh37_OA.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o HipSci_GRCh37_OA.tags_filtered.vcf.gz

#Fix ref and alt alleles
bcftools index HipSci_GRCh37_OA.tags_filtered.vcf.gz
bcftools +fixref HipSci_GRCh37_OA.tags_filtered.vcf.gz -Oz -o HipSci_GRCh37_OA.tags_filtered_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz
mv HipSci_GRCh37_OA.tags_filtered_fixref.vcf.gz HipSci_GRCh37_OA_genotyped.vcf.gz

#Extract chromosomes
bcftools view -r 1 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr1.vcf.gz
bcftools view -r 2 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr2.vcf.gz
bcftools view -r 3 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr3.vcf.gz
bcftools view -r 4 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr4.vcf.gz
bcftools view -r 5 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr5.vcf.gz
bcftools view -r 6 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr6.vcf.gz
bcftools view -r 7 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr7.vcf.gz
bcftools view -r 8 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr8.vcf.gz
bcftools view -r 9 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr9.vcf.gz
bcftools view -r 10 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr10.vcf.gz
bcftools view -r 11 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr11.vcf.gz
bcftools view -r 12 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr12.vcf.gz
bcftools view -r 13 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr13.vcf.gz
bcftools view -r 14 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr14.vcf.gz
bcftools view -r 15 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr15.vcf.gz
bcftools view -r 16 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr16.vcf.gz
bcftools view -r 17 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr17.vcf.gz
bcftools view -r 18 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr18.vcf.gz
bcftools view -r 19 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr19.vcf.gz
bcftools view -r 20 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr20.vcf.gz
bcftools view -r 21 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr21.vcf.gz
bcftools view -r 22 HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chr22.vcf.gz
bcftools view -r X HipSci_GRCh37_OA_genotyped.vcf.gz -Oz -o by_chr/HipSci_GRCh37_chrX.vcf.gz

#Unzip
7za x chr_1.zip -p'6_:TMYXqvs4uFx'
7za x chr_2.zip -p'6_:TMYXqvs4uFx'
7za x chr_3.zip -p'6_:TMYXqvs4uFx'
7za x chr_4.zip -p'6_:TMYXqvs4uFx'
7za x chr_5.zip -p'6_:TMYXqvs4uFx'
7za x chr_6.zip -p'6_:TMYXqvs4uFx'
7za x chr_7.zip -p'6_:TMYXqvs4uFx'
7za x chr_8.zip -p'6_:TMYXqvs4uFx'
7za x chr_9.zip -p'6_:TMYXqvs4uFx'
7za x chr_10.zip -p'6_:TMYXqvs4uFx'
7za x chr_11.zip -p'6_:TMYXqvs4uFx'
7za x chr_12.zip -p'6_:TMYXqvs4uFx'
7za x chr_13.zip -p'6_:TMYXqvs4uFx'
7za x chr_14.zip -p'6_:TMYXqvs4uFx'
7za x chr_15.zip -p'6_:TMYXqvs4uFx'
7za x chr_16.zip -p'6_:TMYXqvs4uFx'
7za x chr_17.zip -p'6_:TMYXqvs4uFx'
7za x chr_18.zip -p'6_:TMYXqvs4uFx'
7za x chr_19.zip -p'6_:TMYXqvs4uFx'
7za x chr_20.zip -p'6_:TMYXqvs4uFx'
7za x chr_21.zip -p'6_:TMYXqvs4uFx'
7za x chr_22.zip -p'6_:TMYXqvs4uFx'


#Controlled access genotypes
plink --file HipSci_test5 --make-bed --out HipSci_GRCh37

#Update build
update_build.sh HipSci_GRCh37 strand_files/HumanCoreExome-12-v1-0-D-b37.strand HipSci_12v1-0_GRCh37

#Convert to VCF
plink --bfile HipSci_12v1-0_GRCh37 --recode vcf-iid --out HipSci_12v1-0_GRCh37

#Identify individuals with high level of missingness
module load vcftools-0.1.13
vcftools --vcf HipSci_12v1-0_GRCh37.vcf --missing-indv --out HipSci_12v1-0_GRCh37_missing