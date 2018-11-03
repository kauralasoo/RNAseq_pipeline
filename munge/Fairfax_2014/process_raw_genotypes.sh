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
update_build.sh Fairfax_2014_set2_XY_sex_merged HumanOmniExpress-12v1_A-b37.strand Fairfax_2014_set2_b37

#Convert to VCF
plink --bfile Fairfax_2014_set2_b37 --recode vcf-iid --out Fairfax_2014_set2
module load bcftools-1.9
bgzip Fairfax_2014_set2.vcf
bcftools index Fairfax_2014_set2.vcf.gz

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
update_build.sh Fairfax_2014_set1_XY_sex_merged HumanOmniExpress-12v1_A-b37.strand Fairfax_2014_set1_b37

#Convert to VCF
plink --bfile Fairfax_2014_set1_b37 --recode vcf-iid --out Fairfax_2014_set1
module load bcftools-1.9
bgzip Fairfax_2014_set1.vcf
bcftools index Fairfax_2014_set1.vcf.gz
