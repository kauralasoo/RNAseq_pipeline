#Convert BED to MAP
plink -bfile ROSMAP_arrayGenotype --recode --tab --out ROSMAP_raw

#Lift over Coordnates to GRCh37
liftOverPlink.py -m ROSMAP_raw.map -p ROSMAP_raw.ped -o ROSMAP_raw_GRCh37 -c NCBI36_to_GRCh37.chain.gz

#Convert back to bed/bim/fam
plink --file ROSMAP_raw_GRCh37 --make-bed --out ROSMAP_raw_GRCh37_bed

#Remove individuals with high levels of missingness
grep -f ~/datasets/controlled_access/SampleArcheology/studies/ROSMAP/ROSMAP_affy_missing_individuals.txt ROSMAP_raw_GRCh37_bed.fam | cut -f1,2 -d " " > remove_list.txt
plink -bfile ROSMAP_raw_GRCh37_bed --remove remove_list.txt --make-bed -out ROSMAP_raw_GRCh37_bed_filtered

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

#Merge both arrays
bcftools merge affy/GRCh38/ROSMAP_affy.vcf.gz illumina/GRCh38/ROSMAP_illumina.vcf.gz -Oz -o ROSMAP_GRCh38_merged.vcf.gz

#Keep RNA-seq individuals
bcftools view -S ~/datasets/controlled_access/SampleArcheology/studies/ROSMAP/ROSMAP_rna_genotype_ids.txt --force-samples ROSMAP_GRCh38_merged.vcf.gz | bcftools filter -i 'F_MISSING < 0.05 & MAF[0] > 0.01' -Oz -o merged/ROSMAP_GRCh38_filtered.vcf.gz 

#Split multi-allelic variants and remove duplicates
bcftools norm -m-any ROSMAP_GRCh38_filtered.vcf.gz | bcftools norm -d all -Oz -o ROSMAP_GRCh38_filtered.no_multi.vcf.gz
