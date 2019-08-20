#Convert BED to MAP
plink -bfile ROSMAP_arrayGenotype --recode --tab --out ROSMAP_raw

#Lift over Coordnates to GRCh37
liftOverPlink.py -m ROSMAP_raw.map -p ROSMAP_raw.ped -o ROSMAP_raw_GRCh37 -c NCBI36_to_GRCh37.chain.gz

#Convert back to bed/bim/fam
plink --file ROSMAP_raw_GRCh37 --make-bed --out ROSMAP_raw_GRCh37_bed

#Harmonize Genotypes
java -jar ~/software/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar --input ROSMAP_raw_GRCh37_bed --inputType PLINK_BED --ref /gpfs/hpc/home/a72094/datasets/1000G/GRCh37_allele_frequencies --refType VCF --update-id --output ROSMAP_harmonized

#Convert PLINK to VCF
module load plink-1.9
plink --bfile ROSMAP_harmonized --recode vcf-iid --out ROSMAP

#Conmpress and index
bgzip ROSMAP.vcf
bcftools index ROSMAP.vcf.gz

#Fix ref and alt alleles
bcftools index ROSMAP.vcf.gz
bcftools +fixref ROSMAP.vcf.gz -Oz -o ROSMAP.fixref.vcf.gz -- -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Remove variants with missing alt alleles
#bcftools filter -e "ALT='.'" ROSMAP.fixref.vcf.gz -Oz -o ROSMAP_GRCh37_sorted.dbSNP.fixref.sorted.noALT.vcf.gz

#Add tags
bcftools +fill-tags ROSMAP.fixref.vcf.gz -Oz -o ROSMAP.fixref.tags.vcf.gz

#Filter rare and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' ROSMAP.fixref.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o ROSMAP.fixref.tags.filtered.vcf.gz

#Remove duplicates and multi-allelics
bcftools norm -d all ROSMAP_tags.filtered.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o ROSMAP_tags.filtered.no_dup.vcf.gz

#Remove non-ref allels
bcftools norm -cx -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa ROSMAP_tags.filtered.no_dup.vcf.gz -Oz -o ROSMAP_GRCh37_genotyped.vcf.gz

#Identify individuals with high level of missing genotypes
module load vcftools-0.1.13
vcftools --gzvcf ROSMAP_GRCh37_genotyped.vcf.gz --missing-indv --out ROSMAP_GRCh37_genotyped

#Extract chromosomes
bcftools index ROSMAP_GRCh37_genotyped.vcf.gz
bcftools view -r 1 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr1.vcf.gz
bcftools view -r 2 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr2.vcf.gz
bcftools view -r 3 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr3.vcf.gz
bcftools view -r 4 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr4.vcf.gz
bcftools view -r 5 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr5.vcf.gz
bcftools view -r 6 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr6.vcf.gz
bcftools view -r 7 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr7.vcf.gz
bcftools view -r 8 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr8.vcf.gz
bcftools view -r 9 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr9.vcf.gz
bcftools view -r 10 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr10.vcf.gz
bcftools view -r 11 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr11.vcf.gz
bcftools view -r 12 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr12.vcf.gz
bcftools view -r 13 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr13.vcf.gz
bcftools view -r 14 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr14.vcf.gz
bcftools view -r 15 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr15.vcf.gz
bcftools view -r 16 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr16.vcf.gz
bcftools view -r 17 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr17.vcf.gz
bcftools view -r 18 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr18.vcf.gz
bcftools view -r 19 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr19.vcf.gz
bcftools view -r 20 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr20.vcf.gz
bcftools view -r 21 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr21.vcf.gz
bcftools view -r 22 ROSMAP_GRCh37_genotyped.vcf.gz -Oz -o by_chr/ROSMAP_GRCh37_chr22.vcf.gz

#Decompress all inputed files
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
bcftools filter -i 'MAF[0] > 0.01' Monocytes_Quach_2016_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Quach_2016_GRCh38.filtered.vcf.gz
bcftools index Quach_2016_GRCh38.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags Quach_2016_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\t%R2\n' | gzip > Quach_2016_GRCh38.variant_information.txt.gz
