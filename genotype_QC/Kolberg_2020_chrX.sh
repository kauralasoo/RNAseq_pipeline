### CEDAR
update_build.sh /gpfs/hpc/projects/Kolberg_2020/genotypes/chrX/CEDAR_raw/CEDAR /gpfs/hpc/projects/Kolberg_2020/genotypes/chrX/CEDAR_raw/HumanOmniExpress-12v1_A-b37.strand CEDAR_b37

### Keep only chromosome 23
plink --bfile CEDAR_b37 --chr 23 --make-bed --out CEDAR_b37_chr23

###### Rename samples and keep only QC-passed samples
#Convert to VCF
plink --bfile CEDAR_b37_chr23 --recode vcf-iid --out CEDAR_b37_chr23

#Change sample names
bcftools reheader -s /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/CEDAR/CEDAR_genotype_name_map.txt CEDAR_b37_chr23.vcf > CEDAR_b37_chr23.renamed.vcf

#Keep QC passed samples
bcftools view CEDAR_b37_chr23.renamed.vcf --samples-file /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/unique_array_samples.txt --force-samples -Oz -o CEDAR_b37_chr23.QC_passed.vcf.gz

#Convert back to plink
plink --vcf CEDAR_b37_chr23.QC_passed.vcf.gz --make-bed --out CEDAR_b37_QC

#Impute sex
plink --bfile CEDAR_b37_QC --impute-sex --make-bed --out CEDAR_b37_QC_sex

#Split PAR region
plink --bfile CEDAR_b37_QC_sex --split-x b37 --make-bed --out CEDAR_b37_QC_split_sex

#Set remaining het haploid genotypes to missing
plink --bfile CEDAR_b37_QC_split_sex --make-bed --set-hh-missing --out CEDAR_b37_QC_no_het

### Keep only chromosome 23
plink --bfile CEDAR_b37_QC_no_het --chr 23 --make-bed --out CEDAR_b37_QC_noPAR

#Rename chromosomes in plink file
sed 's/^23/X/g' CEDAR_b37_QC_noPAR.bim > CEDAR_b37_QC_noPAR.new_bim
mv CEDAR_b37_QC_noPAR.new_bim CEDAR_b37_QC_noPAR.bim

#Run Genotype Harmonizer with 1000G reference
module load java-1.8.0_40
java -jar ~/software/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar\
     --input CEDAR_b37_QC_noPAR\
     --inputType PLINK_BED\
     --ref /gpfs/hpc/projects/genomic_references/1000G/GRCh37/1000G_GRCh37_variant_information\
     --refType VCF\
     --update-id\
     --output CEDAR_b37_QC_noPAR_harmonized

#Extract female samples
plink --bfile CEDAR_b37_QC_noPAR_harmonized --filter-females --make-bed --out CEDAR_females
cut -f1 -d' ' CEDAR_females.fam > CEDAR_females.txt

#Convert to VCF
plink --bfile CEDAR_b37_QC_noPAR_harmonized --recode vcf-iid --out CEDAR_harmonised

#Rename chr23 to X
printf '23\tX\n' > 23_to_x.tsv
bcftools annotate --rename-chrs 23_to_x.tsv CEDAR_harmonised.vcf -o CEDAR_harmonised_chrX.vcf

#Remove missing ALT allele
bcftools filter -e "ALT='.'" CEDAR_harmonised_chrX.vcf -Oz -o CEDAR_harmonised_chrX.noALT.vcf.gz

#Fix the ref allele
bcftools +fixref CEDAR_harmonised_chrX.noALT.vcf.gz -Oz -o CEDAR_harmonised_chrX.noALT.fixref.vcf.gz -- -f /gpfs/hpc/projects/genomic_references/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i /gpfs/hpc/projects/genomic_references/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz





plink --bfile cleaner_fileset --recode vcf-iid --out cleaner_fileset

plink --bfile cleaner_fileset --make-bed --out test
