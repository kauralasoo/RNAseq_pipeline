#Extract genotyped variants from the strand file
cut -f 2,3 humanomni25exome-8v1_a-b37.strand | sort -k1,1 -k2,2n > genotyped_positions.txt

#Extract from vcf
bcftools view -T genotyped_positions.txt Bunt_2015_GRCh37.INFO.vcf.gz -Oz -o Bunt_2015_genotyped.vcf.gz

#Keep variants with INFO == 1
bcftools filter -i "INFO/INFO = '1'" Bunt_2015_genotyped.vcf.gz -Oz -o Bunt_2015_genotyped_INFO1.vcf.gz

#Add all INFO tags
bcftools +fill-tags Bunt_2015_genotyped_INFO1.vcf.gz -Oz -o Bunt_2015_tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Bunt_2015_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o Bunt_2015_tags_filtered.vcf.gz

#Fix ref and alt alleles
bcftools index Bunt_2015_tags_filtered.vcf.gz
bcftools +fixref Bunt_2015_tags_filtered.vcf.gz -Oz -o Bunt_2015_tags_filtered_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa Bunt_2015_tags_filtered_fixref.vcf.gz -Oz -o Bunt_2015_GRCh37_genotyped.vcf.gz

#Sort
bcftools sort Bunt_2015_GRCh37_genotyped.vcf.gz -Oz -o Bunt_2015_GRCh37_genotyped_sorted.vcf.gz
bcftools index Bunt_2015_GRCh37_genotyped_sorted.vcf.gz

#Extract chromosomes
bcftools view -r 1 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr1.vcf.gz
bcftools view -r 2 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr2.vcf.gz
bcftools view -r 3 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr3.vcf.gz
bcftools view -r 4 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr4.vcf.gz
bcftools view -r 5 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr5.vcf.gz
bcftools view -r 6 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr6.vcf.gz
bcftools view -r 7 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr7.vcf.gz
bcftools view -r 8 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr8.vcf.gz
bcftools view -r 9 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr9.vcf.gz
bcftools view -r 10 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr10.vcf.gz
bcftools view -r 11 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr11.vcf.gz
bcftools view -r 12 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr12.vcf.gz
bcftools view -r 13 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr13.vcf.gz
bcftools view -r 14 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr14.vcf.gz
bcftools view -r 15 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr15.vcf.gz
bcftools view -r 16 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr16.vcf.gz
bcftools view -r 17 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr17.vcf.gz
bcftools view -r 18 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr18.vcf.gz
bcftools view -r 19 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr19.vcf.gz
bcftools view -r 20 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr20.vcf.gz
bcftools view -r 21 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr21.vcf.gz
bcftools view -r 22 Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chr22.vcf.gz
bcftools view -r X Bunt_2015_GRCh37_genotyped_sorted.vcf.gz -Oz -o by_chr/Bunt_2015_GRCh37_chrX.vcf.gz

#Extraxt all inputed files
7za x chr_1.zip -p'yB4aqtxFTJVE7J'
7za x chr_2.zip -p'yB4aqtxFTJVE7J'
7za x chr_3.zip -p'yB4aqtxFTJVE7J'
7za x chr_4.zip -p'yB4aqtxFTJVE7J'
7za x chr_5.zip -p'yB4aqtxFTJVE7J'
7za x chr_6.zip -p'yB4aqtxFTJVE7J'
7za x chr_7.zip -p'yB4aqtxFTJVE7J'
7za x chr_8.zip -p'yB4aqtxFTJVE7J'
7za x chr_9.zip -p'yB4aqtxFTJVE7J'
7za x chr_10.zip -p'yB4aqtxFTJVE7J'
7za x chr_11.zip -p'yB4aqtxFTJVE7J'
7za x chr_12.zip -p'yB4aqtxFTJVE7J'
7za x chr_13.zip -p'yB4aqtxFTJVE7J'
7za x chr_14.zip -p'yB4aqtxFTJVE7J'
7za x chr_15.zip -p'yB4aqtxFTJVE7J'
7za x chr_16.zip -p'yB4aqtxFTJVE7J'
7za x chr_17.zip -p'yB4aqtxFTJVE7J'
7za x chr_18.zip -p'yB4aqtxFTJVE7J'
7za x chr_19.zip -p'yB4aqtxFTJVE7J'
7za x chr_20.zip -p'yB4aqtxFTJVE7J'
7za x chr_21.zip -p'yB4aqtxFTJVE7J'
7za x chr_22.zip -p'yB4aqtxFTJVE7J'