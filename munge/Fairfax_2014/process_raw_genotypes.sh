#Extract ped, map and lgen files for set 2
zcat Oxford_eQTL_cohort_genotype2.txt.gz | tail -n+2 | cut -f 2 | uniq | awk '{print $1" "$1" 0 0 0 -9"}' > Fairfax_2014_set2.fam
zcat Oxford_eQTL_cohort_genotype2.txt.gz | tail -n+2 | awk '{print $5" "$1" 0 "$6}' | sort | uniq > Fairfax_2014_set2.map
zcat Oxford_eQTL_cohort_genotype2.txt.gz | tail -n+2 | awk '{print $2" "$2" "$1" "$3" "$4}' | awk '{gsub("-","0",$4); gsub("-","0",$5); print $0}' > Fairfax_2014_set2.lgen

#Convert to binary plink format
plink/plink --lfile Fairfax_2014_set2 --make-bed --out Fairfax_2014_set2_bed

#Imüute sex
plink --bfile Fairfax_2014_set2_bed --impute-sex --make-bed --out Fairfax_2014_set2_sex
plink --bfile Fairfax_2014_set2_sex --merge-x --make-bed --out Fairfax_2014_set2_XY_sex_merged

#Update reference build
update_build.sh Fairfax_2014_set2_XY_sex_merged HumanOmniExpress-12v1_A-b37.Source.strand Fairfax_2014_set2_b37

#Convert to VCF
plink --bfile Fairfax_2014_set2_b37 --recode vcf-iid --out Fairfax_2014_set2
module load bcftools-1.9
bgzip Fairfax_2014_set2.vcf
bcftools index Fairfax_2014_set2.vcf.gz

#Remove missing ALT allele
bcftools filter -e "ALT='.'" Fairfax_2014_set2.vcf.gz -Oz -o Fairfax_2014_set2_noALT.vcf.gz

#Fix the ref allele
bcftools +fixref Fairfax_2014_set2_noALT.vcf.gz -Oz -o Fairfax_2014_set2_noALT_fixref.vcf.gz -- -f ~/hpc/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/hpc/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort the vcf file 
bcftools sort Fairfax_2014_set2_noALT_fixref.vcf.gz -Oz -o Fairfax_2014_set2_noALT_fixref_sorted.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa Fairfax_2014_set2_noALT_fixref_sorted.vcf.gz -Oz -o Fairfax_2014_set2_noALT_fixref_sorted_noref.vcf.gz


#And for set 2
zcat Oxford_eQTL_cohort_genotype.txt.gz | tail -n+2 | cut -f 2 | uniq | awk '{print $1" "$1" 0 0 0 -9"}' > Fairfax_2014_set1.fam
zcat Oxford_eQTL_cohort_genotype.txt.gz | tail -n+2 | awk '{print $5" "$1" 0 "$6}' | sort | uniq > Fairfax_2014_set1.map
zcat Oxford_eQTL_cohort_genotype.txt.gz | tail -n+2 | awk '{print $2" "$2" "$1" "$3" "$4}' | awk '{gsub("-","0",$4); gsub("-","0",$5); print $0}' > Fairfax_2014_set1.lgen

#Convert to binary plink format
plink/plink --lfile Fairfax_2014_set1 --make-bed --out Fairfax_2014_set1_bed

#Imüute sex
plink --bfile Fairfax_2014_set1_bed --impute-sex --make-bed --out Fairfax_2014_set1_sex
plink --bfile Fairfax_2014_set1_sex --merge-x --make-bed --out Fairfax_2014_set1_XY_sex_merged

#Update reference build
#update_build.sh Fairfax_2014_set1_XY_sex_merged HumanOmniExpress-12v1_A-b37.Source.strand Fairfax_2014_set1_b37

#Convert to VCF
plink --bfile Fairfax_2014_set1_XY_sex_merged --recode vcf-iid --out Fairfax_2014_set1
module load bcftools-1.9
bgzip Fairfax_2014_set1.vcf

#Remove missing ALT allele
bcftools filter -e "ALT='.'" Fairfax_2014_set1.vcf.gz -Oz -o Fairfax_2014_set1_noALT.vcf.gz

#Fix the ref allele
bcftools +fixref Fairfax_2014_set1_noALT.vcf.gz -Oz -o Fairfax_2014_set1_noALT_fixref.vcf.gz -- -f ~/hpc/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/hpc/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort the vcf file 
bcftools sort Fairfax_2014_set1_noALT_fixref.vcf.gz -Oz -o Fairfax_2014_set1_noALT_fixref_sorted.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa Fairfax_2014_set1_noALT_fixref_sorted.vcf.gz -Oz -o Fairfax_2014_set1_noALT_fixref_sorted_noref.vcf.gz



#Merge both files
bcftools merge Fairfax_2014_set1_noALT_fixref_sorted_noref.vcf.gz Fairfax_2014_set2_noALT_fixref_sorted_noref.vcf.gz -Oz -o merged.vcf.gz

#Add all INFO tags
bcftools +fill-tags merged.vcf.gz -Oz -o Fairfax_2014_tagged.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Fairfax_2014_tagged.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o Fairfax_2014_filtered.vcf.gz

TODO: Change chr names from 23 and 24 to X and Y

#Remove duplicates and multi-allelics
bcftools norm -d all Fairfax_2014_filtered.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o Fairfax_2014_GRCh37.vcf.gz
bcftools index Fairfax_2014_GRCh37.vcf.gz

#Extract chromosomes
bcftools view -r 1 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr1.vcf.gz
bcftools view -r 2 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr2.vcf.gz
bcftools view -r 3 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr3.vcf.gz
bcftools view -r 4 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr4.vcf.gz
bcftools view -r 5 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr5.vcf.gz
bcftools view -r 6 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr6.vcf.gz
bcftools view -r 7 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr7.vcf.gz
bcftools view -r 8 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr8.vcf.gz
bcftools view -r 9 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr9.vcf.gz
bcftools view -r 10 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr10.vcf.gz
bcftools view -r 11 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr11.vcf.gz
bcftools view -r 12 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr12.vcf.gz
bcftools view -r 13 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr13.vcf.gz
bcftools view -r 14 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr14.vcf.gz
bcftools view -r 15 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr15.vcf.gz
bcftools view -r 16 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr16.vcf.gz
bcftools view -r 17 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr17.vcf.gz
bcftools view -r 18 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr18.vcf.gz
bcftools view -r 19 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr19.vcf.gz
bcftools view -r 20 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr20.vcf.gz
bcftools view -r 21 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr21.vcf.gz
bcftools view -r 22 Fairfax_2014_GRCh37.vcf.gz -Oz -o by_chr/Fairfax_2014_GRCh37_chr22.vcf.gz


