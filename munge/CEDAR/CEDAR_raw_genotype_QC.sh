#Update build
module load plink-1.9.0
update_build.sh CEDAR HumanOmniExpress-12v1_A-b37.strand CEDAR_b37

#Split X
plink --bfile CEDAR_b37 --split-x b37 --make-bed --out CEDAR_b37_XY
plink --bfile CEDAR_b37_XY --impute-sex --make-bed --out CEDAR_b37_XY_imputed
plink --bfile CEDAR_b37_XY_imputed --merge-x --make-bed --out CEDAR_b37_XY_imputed_merged

#Convert to VCF
plink --bfile CEDAR_b37_XY_imputed_merged --recode vcf-iid --output-chr M --out CEDAR_GRCh37

#Count missing variants
module load vcftools-0.1.13
vcftools --vcf CEDAR_GRCh37.vcf --missing-indv --out CEDAR_missing_individuals
vcftools --vcf CEDAR_GRCh37.vcf --missing-site --out CEDAR_missing_sites

#Remove individuals with high missingness (>5% of variants)
bcftools view -S ^CEDAR_exclude_missing.txt CEDAR_GRCh37.vcf -Oz -o CEDAR_GRCh37_nomissing.v^CEDAR_exclude_missing.gz

#Add all INFO tags
bcftools +fill-tags CEDAR_GRCh37_nomissing.vcf.gz -Oz -o CEDAR_GRCh37_nomissing_tags.vgs.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' CEDAR_GRCh37_nomissing_tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o CEDAR_GRCh37_nomissing_tags_filtered.vcf.gz

#Fix the ref allele
bcftools index CEDAR_GRCh37_nomissing_tags_filtered.vcf.gz
bcftools +fixref CEDAR_GRCh37_nomissing_tags_filtered.vcf.gz -Oz -o CEDAR_GRCh37_nomissing_tags_filtered_fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort the vcf file 
bcftools sort CEDAR_GRCh37_nomissing_tags_filtered_fixref.vcf.gz -Oz -o CEDAR_fixref_sorted.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa CEDAR_fixref_sorted.vcf.gz -Oz -o CEDAR_fixref_sorted_noref.vcf.gz

#Check AF distributions
bcftools index CEDAR_fixref_sorted_noref.vcf.gz
bcftools annotate -c INFO/AF DAR_fixref_sorted_noref ~/rocket/datasets/100 -c/GRCh37_allele_frequencies.vcf.gz CEDAR_fixref_sorted_noref.vcf.gz | bcftools +af-dist | grep ^PROB > CEDAR_AFdist.txt

#Remove duplicates and multi-allelics
bcftools norm -d all CEDAR_fixref_sorted_noref.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o CEDAR_fixref_sorted_noref_nodup.vcf.gz

#Rename
mv CEDAR_fixref_sorted_noref_nodup.vcf.gz CEDAR_GRCh37_genotyped.vcf.gz

#Extract chromosomes
bcftools view -r 1 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr1.vcf.gz
bcftools view -r 2 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr2.vcf.gz
bcftools view -r 3 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr3.vcf.gz
bcftools view -r 4 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr4.vcf.gz
bcftools view -r 5 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr5.vcf.gz
bcftools view -r 6 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr6.vcf.gz
bcftools view -r 7 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr7.vcf.gz
bcftools view -r 8 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr8.vcf.gz
bcftools view -r 9 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr9.vcf.gz
bcftools view -r 10 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr10.vcf.gz
bcftools view -r 11 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr11.vcf.gz
bcftools view -r 12 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr12.vcf.gz
bcftools view -r 13 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr13.vcf.gz
bcftools view -r 14 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr14.vcf.gz
bcftools view -r 15 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr15.vcf.gz
bcftools view -r 16 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr16.vcf.gz
bcftools view -r 17 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr17.vcf.gz
bcftools view -r 18 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr18.vcf.gz
bcftools view -r 19 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr19.vcf.gz
bcftools view -r 20 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr20.vcf.gz
bcftools view -r 21 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr21.vcf.gz
bcftools view -r 22 CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chr22.vcf.gz
bcftools view -r X CEDAR_GRCh37_genotyped.vcf.gz -Oz -o by_chr/CEDAR_GRCh37_chrX.vcf.gz


#Extraxt all inputed files
7za x chr_1.zip -p'syh0evL7VUeUAP'
7za x chr_2.zip -p'syh0evL7VUeUAP'
7za x chr_3.zip -p'syh0evL7VUeUAP'
7za x chr_4.zip -p'syh0evL7VUeUAP'
7za x chr_5.zip -p'syh0evL7VUeUAP'
7za x chr_6.zip -p'syh0evL7VUeUAP'
7za x chr_7.zip -p'syh0evL7VUeUAP'
7za x chr_8.zip -p'syh0evL7VUeUAP'
7za x chr_9.zip -p'syh0evL7VUeUAP'
7za x chr_10.zip -p'syh0evL7VUeUAP'
7za x chr_11.zip -p'syh0evL7VUeUAP'
7za x chr_12.zip -p'syh0evL7VUeUAP'
7za x chr_13.zip -p'syh0evL7VUeUAP'
7za x chr_14.zip -p'syh0evL7VUeUAP'
7za x chr_15.zip -p'syh0evL7VUeUAP'
7za x chr_16.zip -p'syh0evL7VUeUAP'
7za x chr_17.zip -p'syh0evL7VUeUAP'
7za x chr_18.zip -p'syh0evL7VUeUAP'
7za x chr_19.zip -p'syh0evL7VUeUAP'
7za x chr_20.zip -p'syh0evL7VUeUAP'
7za x chr_21.zip -p'syh0evL7VUeUAP'
7za x chr_22.zip -p'syh0evL7VUeUAP'

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.01' CEDAR_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o CEDAR_GRCh38.filtered.vcf.gz
bcftools index CEDAR_GRCh38.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags CEDAR_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT\\t%TYPE\\t%AC\\t%AN\\t%MAF\\t%R2\\n' | gzip > CEDAR_GRCh38.variant_information.txt.gz
