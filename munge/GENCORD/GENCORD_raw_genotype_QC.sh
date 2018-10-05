#Index raw genotypes
bcftools index _EGAZ00001014893_GENCORD2.chr10.vcf.gz
bcftools index _EGAZ00001014894_GENCORD2.chr11.vcf.gz
bcftools index _EGAZ00001014895_GENCORD2.chr12.vcf.gz
bcftools index _EGAZ00001014896_GENCORD2.chr13.vcf.gz
bcftools index _EGAZ00001014897_GENCORD2.chr14.vcf.gz
bcftools index _EGAZ00001014898_GENCORD2.chr15.vcf.gz
bcftools index _EGAZ00001014899_GENCORD2.chr16.vcf.gz
bcftools index _EGAZ00001014900_GENCORD2.chr17.vcf.gz
bcftools index _EGAZ00001014901_GENCORD2.chr18.vcf.gz
bcftools index _EGAZ00001014902_GENCORD2.chr19.vcf.gz
bcftools index _EGAZ00001014903_GENCORD2.chr1.vcf.gz
bcftools index _EGAZ00001014904_GENCORD2.chr20.vcf.gz
bcftools index _EGAZ00001014905_GENCORD2.chr21.vcf.gz
bcftools index _EGAZ00001014906_GENCORD2.chr22.vcf.gz
bcftools index _EGAZ00001014907_GENCORD2.chr2.vcf.gz
bcftools index _EGAZ00001014908_GENCORD2.chr3.vcf.gz
bcftools index _EGAZ00001014909_GENCORD2.chr4.vcf.gz
bcftools index _EGAZ00001014910_GENCORD2.chr5.vcf.gz
bcftools index _EGAZ00001014911_GENCORD2.chr6.vcf.gz
bcftools index _EGAZ00001014912_GENCORD2.chr7.vcf.gz
bcftools index _EGAZ00001014913_GENCORD2.chr8.vcf.gz
bcftools index _EGAZ00001014914_GENCORD2.chr9.vcf.gz
bcftools index _EGAZ00001014915_GENCORD2.chrX.vcf.gz

#Concat raw files
bcftools concat _EGAZ00001014903_GENCORD2.chr1.vcf.gz _EGAZ00001014907_GENCORD2.chr2.vcf.gz _EGAZ00001014908_GENCORD2.chr3.vcf.gz _EGAZ00001014909_GENCORD2.chr4.vcf.gz _EGAZ00001014910_GENCORD2.chr5.vcf.gz _EGAZ00001014911_GENCORD2.chr6.vcf.gz _EGAZ00001014912_GENCORD2.chr7.vcf.gz _EGAZ00001014913_GENCORD2.chr8.vcf.gz _EGAZ00001014914_GENCORD2.chr9.vcf.gz _EGAZ00001014893_GENCORD2.chr10.vcf.gz _EGAZ00001014894_GENCORD2.chr11.vcf.gz _EGAZ00001014895_GENCORD2.chr12.vcf.gz _EGAZ00001014896_GENCORD2.chr13.vcf.gz _EGAZ00001014897_GENCORD2.chr14.vcf.gz _EGAZ00001014898_GENCORD2.chr15.vcf.gz _EGAZ00001014899_GENCORD2.chr16.vcf.gz _EGAZ00001014900_GENCORD2.chr17.vcf.gz _EGAZ00001014901_GENCORD2.chr18.vcf.gz _EGAZ00001014902_GENCORD2.chr19.vcf.gz _EGAZ00001014904_GENCORD2.chr20.vcf.gz _EGAZ00001014905_GENCORD2.chr21.vcf.gz _EGAZ00001014906_GENCORD2.chr22.vcf.gz _EGAZ00001014915_GENCORD2.chrX.vcf.gz -Oz -o GENCORD_raw_genotypes.vcf.gz

#Keep only variants that were directly genotyped
bcftools filter -i 'INFO/GT = "1"' GENCORD_raw_genotypes.vcf.gz -Oz -o GENCORD_genotyped.vcf.gz

#Set imputed variants to missing
bcftools +setGT GENCORD_genotyped.vcf.gz -Oz -o GENCORD_genotyped.set_missing.vcf.gz -- --target-gt 'q' --new-gt '.' --include 'GQ < 1'

#Remove GQ field
bcftools annotate -x FORMAT/GQ GENCORD_genotyped.set_missing.vcf.gz -Oz -o GENCORD_genotyped.set_missing.noGQ.vcf.gz

#Identify individuals with high level of missing genotypes
module load vcftools-0.1.13s
vcftools --gzvcf GENCORD_genotyped.set_missing.noGQ.vcf.gz --missing-indv --out GENCORD_missing

#Convert to plink and back
plink --vcf GENCORD_genotyped.set_missing.noGQ.vcf.gz --make-bed --out GENCORD
plink --bfile GENCORD --recode vcf-iid --out GENCORD_genotyped.plinked

#Add all INFO tags
bcftools +fill-tags GENCORD_genotyped.plinked.vcf -Oz -o GENCORD_tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' GENCORD_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o GENCORD_tags_filtered.vcf.gz

#Remane chromosomes
bcftools annotate --rename-chrs ../chromsome_map.txt GENCORD_tags_filtered.vcf.gz -Oz -o GENCORD_tags_filtered_rename.vcf.gz

#Fix ref and alt alleles
bcftools index GENCORD_tags_filtered_rename.vcf.gz
bcftools +fixref GENCORD_tags_filtered_rename.vcf.gz -Oz -o GENCORD_tags_filtered_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa GENCORD_tags_filtered_fixref.vcf.gz -Oz -o GENCORD_GRCh37_genotyped.vcf.gz
bcftools index GENCORD_GRCh37_genotyped.vcf.gz

#Extract chromosomes
bcftools view -r 1 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr1.vcf.gz
bcftools view -r 2 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr2.vcf.gz
bcftools view -r 3 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr3.vcf.gz
bcftools view -r 4 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr4.vcf.gz
bcftools view -r 5 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr5.vcf.gz
bcftools view -r 6 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr6.vcf.gz
bcftools view -r 7 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr7.vcf.gz
bcftools view -r 8 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr8.vcf.gz
bcftools view -r 9 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr9.vcf.gz
bcftools view -r 10 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr10.vcf.gz
bcftools view -r 11 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr11.vcf.gz
bcftools view -r 12 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr12.vcf.gz
bcftools view -r 13 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr13.vcf.gz
bcftools view -r 14 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr14.vcf.gz
bcftools view -r 15 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr15.vcf.gz
bcftools view -r 16 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr16.vcf.gz
bcftools view -r 17 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr17.vcf.gz
bcftools view -r 18 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr18.vcf.gz
bcftools view -r 19 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr19.vcf.gz
bcftools view -r 20 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr20.vcf.gz
bcftools view -r 21 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr21.vcf.gz
bcftools view -r 22 GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chr22.vcf.gz
bcftools view -r X GENCORD_GRCh37_genotyped.vcf.gz -Oz -o by_chr/GENCORD_GRCh37_chrX.vcf.gz

#Extraxt all inputed files
7za x chr_1.zip -p'Thr6tAkSFIQfU2'
7za x chr_2.zip -p'Thr6tAkSFIQfU2'
7za x chr_3.zip -p'Thr6tAkSFIQfU2'
7za x chr_4.zip -p'Thr6tAkSFIQfU2'
7za x chr_5.zip -p'Thr6tAkSFIQfU2'
7za x chr_6.zip -p'Thr6tAkSFIQfU2'
7za x chr_7.zip -p'Thr6tAkSFIQfU2'
7za x chr_8.zip -p'Thr6tAkSFIQfU2'
7za x chr_9.zip -p'Thr6tAkSFIQfU2'
7za x chr_10.zip -p'Thr6tAkSFIQfU2'
7za x chr_11.zip -p'Thr6tAkSFIQfU2'
7za x chr_12.zip -p'Thr6tAkSFIQfU2'
7za x chr_13.zip -p'Thr6tAkSFIQfU2'
7za x chr_14.zip -p'Thr6tAkSFIQfU2'
7za x chr_15.zip -p'Thr6tAkSFIQfU2'
7za x chr_16.zip -p'Thr6tAkSFIQfU2'
7za x chr_17.zip -p'Thr6tAkSFIQfU2'
7za x chr_18.zip -p'Thr6tAkSFIQfU2'
7za x chr_19.zip -p'Thr6tAkSFIQfU2'
7za x chr_20.zip -p'Thr6tAkSFIQfU2'
7za x chr_21.zip -p'Thr6tAkSFIQfU2'
7za x chr_22.zip -p'Thr6tAkSFIQfU2'

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.01' GENCORD_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o GENCORD_GRCh38.filtered.vcf.gz
bcftools index GENCORD_GRCh38.filtered.vcf.gz
