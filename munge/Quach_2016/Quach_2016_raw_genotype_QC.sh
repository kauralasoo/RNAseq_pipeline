#Convert PLINK to VCF
module load plink-1.9
plink --bfile EvoImmunoPop_Omni5exome_200x3773685 --recode vcf-iid --out Quach_2016

#Add SNP ids from dbSNP
bgzip Quach_2016.vcf
bcftools index Quach_2016.vcf.gz
bcftools annotate -a /gpfs/hpchome/a72094/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz -c ID Quach_2016.vcf.gz -Oz -o Quach_2016.dbSNP.vcf.gz

#Fix ref and alt alleles
bcftools index Quach_2016.dbSNP.vcf.gz
bcftools +fixref Quach_2016.dbSNP.vcf.gz -Oz -o Quach_2016.dbSNP.fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Remove variants with missing alt alleles
bcftools filter -e "ALT='.'" Quach_2016.dbSNP.fixref.vcf.gz -Oz -o Quach_2016.dbSNP.fixref.noALT.vcf.gz

#Add tags
bcftools +fill-tags Quach_2016.dbSNP.fixref.noALT.vcf.gz -Oz -o Quach_2016.dbSNP.fixref.noALT.tags.vcf.gz

#Filter rare and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Quach_2016.dbSNP.fixref.noALT.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o Quach_2016.dbSNP.fixref.noALT.tags.filtered.vcf.gz

#Sort the vcf file 
bcftools sort Quach_2016.dbSNP.fixref.noALT.tags.filtered.vcf.gz -Oz -o Quach_2016.dbSNP.fixref.noALT.tags.filtered.sorted.vcf.gz

#Remove duplicates and multi-allelics
bcftools norm -d all Quach_2016.dbSNP.fixref.noALT.tags.filtered.sorted.vcf.gz | bcftools norm -m+any | bcftools view -m2 -M2 -Oz -o Quach_2016.no_dup.vcf.gz

#Remove non-ref allels
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa Quach_2016.no_dup.vcf.gz -Oz -o Quach_2016_GRCh37_genotyped.vcf.gz

#Identify individuals with high level of missing genotypes
module load vcftools-0.1.13
vcftools --gzvcf Quach_2016_GRCh37_genotyped.vcf.gz --missing-indv --out Quach_2016_GRCh37_genotyped

#Extract chromosomes
bcftools index Quach_2016_GRCh37_genotyped.vcf.gz
bcftools view -r 1 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr1.vcf.gz
bcftools view -r 2 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr2.vcf.gz
bcftools view -r 3 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr3.vcf.gz
bcftools view -r 4 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr4.vcf.gz
bcftools view -r 5 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr5.vcf.gz
bcftools view -r 6 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr6.vcf.gz
bcftools view -r 7 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr7.vcf.gz
bcftools view -r 8 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr8.vcf.gz
bcftools view -r 9 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr9.vcf.gz
bcftools view -r 10 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr10.vcf.gz
bcftools view -r 11 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr11.vcf.gz
bcftools view -r 12 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr12.vcf.gz
bcftools view -r 13 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr13.vcf.gz
bcftools view -r 14 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr14.vcf.gz
bcftools view -r 15 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr15.vcf.gz
bcftools view -r 16 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr16.vcf.gz
bcftools view -r 17 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr17.vcf.gz
bcftools view -r 18 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr18.vcf.gz
bcftools view -r 19 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr19.vcf.gz
bcftools view -r 20 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr20.vcf.gz
bcftools view -r 21 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr21.vcf.gz
bcftools view -r 22 Quach_2016_GRCh37_genotyped.vcf.gz -Oz -o by_chr/Quach_2016_GRCh37_chr22.vcf.gz

#Decompress all inputed files
7za x chr_1.zip -p'OcBXHL05bRnvs'
7za x chr_2.zip -p'OcBXHL05bRnvs'
7za x chr_3.zip -p'OcBXHL05bRnvs'
7za x chr_4.zip -p'OcBXHL05bRnvs'
7za x chr_5.zip -p'OcBXHL05bRnvs'
7za x chr_6.zip -p'OcBXHL05bRnvs'
7za x chr_7.zip -p'OcBXHL05bRnvs'
7za x chr_8.zip -p'OcBXHL05bRnvs'
7za x chr_9.zip -p'OcBXHL05bRnvs'
7za x chr_10.zip -p'OcBXHL05bRnvs'
7za x chr_11.zip -p'OcBXHL05bRnvs'
7za x chr_12.zip -p'OcBXHL05bRnvs'
7za x chr_13.zip -p'OcBXHL05bRnvs'
7za x chr_14.zip -p'OcBXHL05bRnvs'
7za x chr_15.zip -p'OcBXHL05bRnvs'
7za x chr_16.zip -p'OcBXHL05bRnvs'
7za x chr_17.zip -p'OcBXHL05bRnvs'
7za x chr_18.zip -p'OcBXHL05bRnvs'
7za x chr_19.zip -p'OcBXHL05bRnvs'
7za x chr_20.zip -p'OcBXHL05bRnvs'
7za x chr_21.zip -p'OcBXHL05bRnvs'
7za x chr_22.zip -p'OcBXHL05bRnvs'

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.05' Monocytes_Quach_2016_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Quach_2016_GRCh38.filtered.vcf.gz
bcftools index Quach_2016_GRCh38.filtered.vcf.gz