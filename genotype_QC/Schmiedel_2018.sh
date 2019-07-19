#Remove chromosome 0
plink --bfile IC_DNA --not-chr 0 --make-bed --out DICE_no_zero

#Impute sex
plink --bfile DICE_no_zero --impute-sex --make-bed --out DICE_imputed

#Merge X
plink --bfile DICE_imputed --merge-x --make-bed --out DICE_imputed_merged

#Remove duplicate variants
cut -f 2 DICE_imputed_merged.bim | sort | uniq -d > 1.dups
plink --bfile DICE_imputed_merged --exclude 1.dups --make-bed --out DICE_imputed_merged_nodups

#Udate build
update_build.sh DICE_imputed_merged_nodups Multi-EthnicGlobal_D1-b37.strand DICE_build

#Convert to VCF
plink --bfile DICE_build --recode vcf-iid --output-chr M --out DICE_GRCh37

#Count missing variants
module load vcftools-0.1.13
vcftools --vcf DICE_GRCh37.vcf --missing-indv --out DICE_missing_individuals
vcftools --vcf DICE_GRCh37.vcf --missing-site --out DICE_missing_sites

#Add all INFO tags
bcftools +fill-tags DICE_GRCh37.vcf -Oz -o DICE_GRCh37_tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' DICE_GRCh37_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o DICE_GRCh37_tags_filtered.vcf.gz

#Fix the ref allele
bcftools +fixref DICE_GRCh37_tags_filtered.vcf.gz -Oz -o DICE_GRCh37_tags_filtered_fixref.vcf.gz -- -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort the vcf file 
bcftools sort DICE_GRCh37_tags_filtered_fixref.vcf.gz -Oz -o DICE_GRCh37_tags_filtered_fixref_sorted.vcf.gz
bcftools index DICE_GRCh37_tags_filtered_fixref_sorted.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa DICE_GRCh37_tags_filtered_fixref_sorted.vcf.gz -Oz -o DICE_GRCh37_tags_filtered_fixref_sorted_noref.vcf.gz

#Check AF distributions
bcftools index CEDAR_fixref_sorted_noref.vcf.gz
bcftools annotate -c INFO/AF DAR_fixref_sorted_noref ~/rocket/datasets/100 -c/GRCh37_allele_frequencies.vcf.gz CEDAR_fixref_sorted_noref.vcf.gz | bcftools +af-dist | grep ^PROB > CEDAR_AFdist.txt

#Remove duplicates and multi-allelics
bcftools norm -d all DICE_GRCh37_tags_filtered_fixref_sorted.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o DICE_GRCh37_tags_filtered_fixref_sorted_nodup.vcf.gz

#Rename
mv DICE_GRCh37_tags_filtered_fixref_sorted_nodup.vcf.gz DICE_GRCh37_genotyped.vcf.gz
bcftools index DICE_GRCh37_genotyped.vcf.gz

#Extract chromosomes
bcftools view -r 1 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr1.vcf.gz
bcftools view -r 2 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr2.vcf.gz
bcftools view -r 3 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr3.vcf.gz
bcftools view -r 4 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr4.vcf.gz
bcftools view -r 5 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr5.vcf.gz
bcftools view -r 6 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr6.vcf.gz
bcftools view -r 7 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr7.vcf.gz
bcftools view -r 8 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr8.vcf.gz
bcftools view -r 9 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr9.vcf.gz
bcftools view -r 10 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr10.vcf.gz
bcftools view -r 11 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr11.vcf.gz
bcftools view -r 12 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr12.vcf.gz
bcftools view -r 13 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr13.vcf.gz
bcftools view -r 14 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr14.vcf.gz
bcftools view -r 15 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr15.vcf.gz
bcftools view -r 16 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr16.vcf.gz
bcftools view -r 17 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr17.vcf.gz
bcftools view -r 18 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr18.vcf.gz
bcftools view -r 19 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr19.vcf.gz
bcftools view -r 20 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr20.vcf.gz
bcftools view -r 21 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr21.vcf.gz
bcftools view -r 22 DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chr22.vcf.gz
bcftools view -r X DICE_GRCh37_genotyped.vcf.gz -Oz -o by_chr/DICE_GRCh37_chrX.vcf.gz


#Extraxt all inputed files
7za x chr_1.zip -p'password'
7za x chr_2.zip -p'password'
7za x chr_3.zip -p'password'
7za x chr_4.zip -p'password'
7za x chr_5.zip -p'password'
7za x chr_6.zip -p'password'
7za x chr_7.zip -p'password'
7za x chr_8.zip -p'password'
7za x chr_9.zip -p'password'
7za x chr_10.zip -p'password'
7za x chr_11.zip -p'password'
7za x chr_12.zip -p'password'
7za x chr_13.zip -p'password'
7za x chr_14.zip -p'password'
7za x chr_15.zip -p'password'
7za x chr_16.zip -p'password'
7za x chr_17.zip -p'password'
7za x chr_18.zip -p'password'
7za x chr_19.zip -p'password'
7za x chr_20.zip -p'password'
7za x chr_21.zip -p'password'
7za x chr_22.zip -p'password'

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.01' Schmiedel_2018_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Schmiedel_2018_GRCh38.filtered.vcf.gz
bcftools index Schmiedel_2018_GRCh38.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags Schmiedel_2018_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\t%R2\n' | gzip > Schmiedel_2018_GRCh38.variant_information.txt.gz
