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
module load plink-1.9
plink --file controlled/HipSci_v1-0_PLINK/HipSci_v1-0 --make-bed --out controlled/HipSci_v1-0_PLINK/HipSci_GRCh37_v1-0
plink --file controlled/HipSci_v1-1_PLINK/HipSci_v1-1 --make-bed --out controlled/HipSci_v1-1_PLINK/HipSci_GRCh37_v1-1

#Update build
update_build.sh controlled/HipSci_v1-0_PLINK/HipSci_GRCh37_v1-0 controlled/strand_files/HumanCoreExome-12-v1-0-D-b37.strand controlled/processed/HipSci_GRCh37_v1-0
update_build.sh controlled/HipSci_v1-1_PLINK/HipSci_GRCh37_v1-1 controlled/strand_files/HumanCoreExome-12-v1-1-C-b37.strand controlled/processed/HipSci_GRCh37_v1-1

#Convert to VCF
plink --bfile controlled/processed/HipSci_GRCh37_v1-0 --recode vcf-iid --out controlled/processed/HipSci_GRCh37_v1-0
plink --bfile controlled/processed/HipSci_GRCh37_v1-1 --recode vcf-iid --out controlled/processed/HipSci_GRCh37_v1-1

#Identify individuals with high level of missingness
module load vcftools-0.1.13
vcftools --vcf controlled/processed/HipSci_GRCh37_v1-0.vcf --missing-indv --out controlled/processed/HipSci_GRCh37_v1-0
vcftools --vcf controlled/processed/HipSci_GRCh37_v1-1.vcf --missing-indv --out controlled/processed/HipSci_GRCh37_v1-1

#Remove samples with high level of missingness (>5%)
bcftools view controlled/processed/HipSci_GRCh37_v1-0.vcf -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.vcf.gz
bcftools view -s ^HPSI0513pf-dovq,HPSI0713pf-uimo,HPSI0713i-uimo_2 controlled/processed/HipSci_GRCh37_v1-1.vcf -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.vcf.gz

#Add all INFO tags
bcftools +fill-tags controlled/processed/HipSci_GRCh37_v1-0.imiss.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.vcf.gz
bcftools +fill-tags controlled/processed/HipSci_GRCh37_v1-1.imiss.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.vcf.gz
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.vcf.gz

#Remane chromosomes and remove MT sequences
bcftools annotate --rename-chrs ~/rocket/projects/GEUVADIS/RNAseq_pipeline/genotype_scripts/chromsome_names.txt controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.vcf.gz | bcftools view -t ^MT -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.vcf.gz
bcftools annotate --rename-chrs ~/rocket/projects/GEUVADIS/RNAseq_pipeline/genotype_scripts/chromsome_names.txt controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.vcf.gz | bcftools view -t ^MT -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.vcf.gz

#Index
bcftools index controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.vcf.gz
bcftools index controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.vcf.gz

#Add SNP ids from dbSNP
bcftools annotate -a /gpfs/hpchome/a72094/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz -c ID controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.vcf.gz
bcftools annotate -a /gpfs/hpchome/a72094/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz -c ID controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.vcf.gz

#Index
bcftools index controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.vcf.gz
bcftools index controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.vcf.gz

#Fix ref and alt alleles
bcftools +fixref controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz
bcftools +fixref controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.fixref.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

#Sort
bcftools sort controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.fixref.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.fixref.sorted.vcf.gz
bcftools sort controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.fixref.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.fixref.sorted.vcf.gz

#Index
bcftools index controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.fixref.sorted.vcf.gz
bcftools index controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.fixref.sorted.vcf.gz

#Remove remaining non-ref alleles
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa controlled/processed/HipSci_GRCh37_v1-0.imiss.tags.filtered.chr.dbSNP.fixref.sorted.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-0_genotyped.vcf.gz
bcftools norm -cx -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa controlled/processed/HipSci_GRCh37_v1-1.imiss.tags.filtered.chr.dbSNP.fixref.sorted.vcf.gz -Oz -o controlled/processed/HipSci_GRCh37_v1-1_genotyped.vcf.gz

HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz

#Extract chromosomes
bcftools view -r 1 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr1.vcf.gz
bcftools view -r 2 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr2.vcf.gz
bcftools view -r 3 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr3.vcf.gz
bcftools view -r 4 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr4.vcf.gz
bcftools view -r 5 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr5.vcf.gz
bcftools view -r 6 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr6.vcf.gz
bcftools view -r 7 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr7.vcf.gz
bcftools view -r 8 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr8.vcf.gz
bcftools view -r 9 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr9.vcf.gz
bcftools view -r 10 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr10.vcf.gz
bcftools view -r 11 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr11.vcf.gz
bcftools view -r 12 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr12.vcf.gz
bcftools view -r 13 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr13.vcf.gz
bcftools view -r 14 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr14.vcf.gz
bcftools view -r 15 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr15.vcf.gz
bcftools view -r 16 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr16.vcf.gz
bcftools view -r 17 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr17.vcf.gz
bcftools view -r 18 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr18.vcf.gz
bcftools view -r 19 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr19.vcf.gz
bcftools view -r 20 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr20.vcf.gz
bcftools view -r 21 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr21.vcf.gz
bcftools view -r 22 HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chr22.vcf.gz
bcftools view -r X HipSci_GRCh37_v1-1/HipSci_GRCh37_v1-1_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-1/by_chr/HipSci_GRCh37_v1-1_chrX.vcf.gz

bcftools view -r 1 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr1.vcf.gz
bcftools view -r 2 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr2.vcf.gz
bcftools view -r 3 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr3.vcf.gz
bcftools view -r 4 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr4.vcf.gz
bcftools view -r 5 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr5.vcf.gz
bcftools view -r 6 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr6.vcf.gz
bcftools view -r 7 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr7.vcf.gz
bcftools view -r 8 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr8.vcf.gz
bcftools view -r 9 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr9.vcf.gz
bcftools view -r 10 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr10.vcf.gz
bcftools view -r 11 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr11.vcf.gz
bcftools view -r 12 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr12.vcf.gz
bcftools view -r 13 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr13.vcf.gz
bcftools view -r 14 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr14.vcf.gz
bcftools view -r 15 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr15.vcf.gz
bcftools view -r 16 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr16.vcf.gz
bcftools view -r 17 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr17.vcf.gz
bcftools view -r 18 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr18.vcf.gz
bcftools view -r 19 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr19.vcf.gz
bcftools view -r 20 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr20.vcf.gz
bcftools view -r 21 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr21.vcf.gz
bcftools view -r 22 HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chr22.vcf.gz
bcftools view -r X HipSci_GRCh37_v1-0/HipSci_GRCh37_v1-0_genotyped.vcf.gz -Oz -o HipSci_GRCh37_v1-0/by_chr/HipSci_GRCh37_v1-0_chrX.vcf.gz

#Unzip
7za x chr_1.zip -p'Ay(;sn2BAmgMU4'
7za x chr_2.zip -p'Ay(;sn2BAmgMU4'
7za x chr_3.zip -p'Ay(;sn2BAmgMU4'
7za x chr_4.zip -p'Ay(;sn2BAmgMU4'
7za x chr_5.zip -p'Ay(;sn2BAmgMU4'
7za x chr_6.zip -p'Ay(;sn2BAmgMU4'
7za x chr_7.zip -p'Ay(;sn2BAmgMU4'
7za x chr_8.zip -p'Ay(;sn2BAmgMU4'
7za x chr_9.zip -p'Ay(;sn2BAmgMU4'
7za x chr_10.zip -p'Ay(;sn2BAmgMU4'
7za x chr_11.zip -p'Ay(;sn2BAmgMU4'
7za x chr_12.zip -p'Ay(;sn2BAmgMU4'
7za x chr_13.zip -p'Ay(;sn2BAmgMU4'
7za x chr_14.zip -p'Ay(;sn2BAmgMU4'
7za x chr_15.zip -p'Ay(;sn2BAmgMU4'
7za x chr_16.zip -p'Ay(;sn2BAmgMU4'
7za x chr_17.zip -p'Ay(;sn2BAmgMU4'
7za x chr_18.zip -p'Ay(;sn2BAmgMU4'
7za x chr_19.zip -p'Ay(;sn2BAmgMU4'
7za x chr_20.zip -p'Ay(;sn2BAmgMU4'
7za x chr_21.zip -p'Ay(;sn2BAmgMU4'
7za x chr_22.zip -p'Ay(;sn2BAmgMU4'

7za x chr_1.zip -p'2RgCcPHx9sOu8)'
7za x chr_2.zip -p'2RgCcPHx9sOu8)'
7za x chr_3.zip -p'2RgCcPHx9sOu8)'
7za x chr_4.zip -p'2RgCcPHx9sOu8)'
7za x chr_5.zip -p'2RgCcPHx9sOu8)'
7za x chr_6.zip -p'2RgCcPHx9sOu8)'
7za x chr_7.zip -p'2RgCcPHx9sOu8)'
7za x chr_8.zip -p'2RgCcPHx9sOu8)'
7za x chr_9.zip -p'2RgCcPHx9sOu8)'
7za x chr_10.zip -p'2RgCcPHx9sOu8)'
7za x chr_11.zip -p'2RgCcPHx9sOu8)'
7za x chr_12.zip -p'2RgCcPHx9sOu8)'
7za x chr_13.zip -p'2RgCcPHx9sOu8)'
7za x chr_14.zip -p'2RgCcPHx9sOu8)'
7za x chr_15.zip -p'2RgCcPHx9sOu8)'
7za x chr_16.zip -p'2RgCcPHx9sOu8)'
7za x chr_17.zip -p'2RgCcPHx9sOu8)'
7za x chr_18.zip -p'2RgCcPHx9sOu8)'
7za x chr_19.zip -p'2RgCcPHx9sOu8)'
7za x chr_20.zip -p'2RgCcPHx9sOu8)'
7za x chr_21.zip -p'2RgCcPHx9sOu8)'
7za x chr_22.zip -p'2RgCcPHx9sOu8)'

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.01' HipSci_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o HipSci_GRCh38.filtered.vcf.gz
bcftools index HipSci_GRCh38.filtered.vcf.gz