### Pre-imputation QC

### Kasela_2017
#Convert to binary plink format
plink --file ~/datasets/controlled_access/Kasela_2017/genotypes/EGV_raw/CTG_omniSNP_unik313ind --make-bed --out Kasela_2017/Kasela_2017_raw

#Run pre-imputation QC
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile Kasela_2017/Kasela_2017_raw\
 --output_name Kasela_2017_GRCh37_genotyped\
 --outdir preimpute

 ### CEDAR
update_build.sh /gpfs/hpc/home/a72094/datasets/open_access/CEDAR/genotypes/PLINK_100718_1018/CEDAR /gpfs/hpc/home/a72094/datasets/open_access/CEDAR/genotypes/PLINK_100718_1018/HumanOmniExpress-12v1_A-b37.strand CEDAR_b37
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile CEDAR_b37\
 --output_name CEDAR_GRCh37_genotyped\
 --outdir preimpute

#Remove ref mismatch alleles
bcftools norm CEDAR_GRCh37_genotyped.vcf.gz --check-ref x -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -Oz -o CEDAR_GRCh37_genotyped.filtered.vcf.gz

#Change sample names
bcftools reheader -s /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/CEDAR/CEDAR_genotype_name_map.txt CEDAR_GRCh37_genotyped.filtered.vcf.gz > CEDAR_GRCh37_genotyped.renamed.vcf.gz

### Fairfax_2014
#Update reference build
nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/Fairfax_2014/genotypes/processed/Fairfax_2014_set1_b37\
 --output_name Fairfax_2014_set1_GRCh37_genotyped\
 --outdir preimpute

 nextflow run pre-imputation_qc.nf -profile eqtl_catalogue -resume\
 --bfile /gpfs/hpc/home/a72094/datasets/controlled_access/Fairfax_2014/genotypes/plink_raw/Fairfax_2014_set2_bed\
 --output_name Fairfax_2014_set2_GRCh37_genotyped\
 --outdir preimpute

bcftools reheader -s /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/studies/Fairfax_2014/genotype_name_map.txt Fairfax_2014_set1_GRCh37_genotyped.vcf.gz > Fairfax_2014_set1_GRCh37_genotyped.renamed.vcf.gz
bcftools reheader -s /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/studies/Fairfax_2014/genotype_name_map.txt Fairfax_2014_set2_GRCh37_genotyped.vcf.gz > Fairfax_2014_set2_GRCh37_genotyped.renamed.vcf.gz

#Remove non-ref alleles
bcftools norm Kasela_2017_GRCh37_genotyped.vcf.gz --check-ref x -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -Oz -o Kasela_2017_GRCh37_genotyped.filtered.vcf.gz
bcftools norm Fairfax_2014_set1_GRCh37_genotyped.renamed.vcf.gz --check-ref x -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -Oz -o Fairfax_2014_set1_GRCh37_genotyped.filtered.vcf.gz
bcftools norm Fairfax_2014_set2_GRCh37_genotyped.renamed.vcf.gz --check-ref x -f ~/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -Oz -o Fairfax_2014_set2_GRCh37_genotyped.filtered.vcf.gz

#Merge all of the VCF files together
bcftools index CEDAR_GRCh37_genotyped.renamed.vcf.gz
bcftools index Kasela_2017_GRCh37_genotyped.filtered.vcf.gz
bcftools index Fairfax_2014_set2_GRCh37_genotyped.filtered.vcf.gz
bcftools index Fairfax_2014_set1_GRCh37_genotyped.filtered.vcf.gz
bcftools merge CEDAR_GRCh37_genotyped.renamed.vcf.gz Kasela_2017_GRCh37_genotyped.filtered.vcf.gz Fairfax_2014_set1_GRCh37_genotyped.filtered.vcf.gz Fairfax_2014_set2_GRCh37_genotyped.filtered.vcf.gz -Oz -o merged.vcf.gz

#Keep QCd samples
bcftools view merged.vcf.gz --samples-file ~/datasets/controlled_access/SampleArcheology/studies/unique_array_samples.txt -Oz -o filtered.vcf.gz

#Fill tags
bcftools +fill-tags filtered.vcf.gz -Oz -o tagged.vcf.gz

#Final filtering
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' tagged.vcf.gz |\
     bcftools norm -d all |\
     bcftools norm -m+any |\
     bcftools view -m2 -M2 -Oz -o final_filtered.vcf.gz

#Extract chromosomes
bcftools view -r 1 final_filtered.vcf.gz -Oz -o by_chr/chr_1.vcf.gz
bcftools view -r 2 final_filtered.vcf.gz -Oz -o by_chr/chr_2.vcf.gz
bcftools view -r 3 final_filtered.vcf.gz -Oz -o by_chr/chr_3.vcf.gz
bcftools view -r 4 final_filtered.vcf.gz -Oz -o by_chr/chr_4.vcf.gz
bcftools view -r 5 final_filtered.vcf.gz -Oz -o by_chr/chr_5.vcf.gz
bcftools view -r 6 final_filtered.vcf.gz -Oz -o by_chr/chr_6.vcf.gz
bcftools view -r 7 final_filtered.vcf.gz -Oz -o by_chr/chr_7.vcf.gz
bcftools view -r 8 final_filtered.vcf.gz -Oz -o by_chr/chr_8.vcf.gz
bcftools view -r 9 final_filtered.vcf.gz -Oz -o by_chr/chr_9.vcf.gz
bcftools view -r 10 final_filtered.vcf.gz -Oz -o by_chr/chr_10.vcf.gz
bcftools view -r 11 final_filtered.vcf.gz -Oz -o by_chr/chr_11.vcf.gz
bcftools view -r 12 final_filtered.vcf.gz -Oz -o by_chr/chr_12.vcf.gz
bcftools view -r 13 final_filtered.vcf.gz -Oz -o by_chr/chr_13.vcf.gz
bcftools view -r 14 final_filtered.vcf.gz -Oz -o by_chr/chr_14.vcf.gz
bcftools view -r 15 final_filtered.vcf.gz -Oz -o by_chr/chr_15.vcf.gz
bcftools view -r 16 final_filtered.vcf.gz -Oz -o by_chr/chr_16.vcf.gz
bcftools view -r 17 final_filtered.vcf.gz -Oz -o by_chr/chr_17.vcf.gz
bcftools view -r 18 final_filtered.vcf.gz -Oz -o by_chr/chr_18.vcf.gz
bcftools view -r 19 final_filtered.vcf.gz -Oz -o by_chr/chr_19.vcf.gz
bcftools view -r 20 final_filtered.vcf.gz -Oz -o by_chr/chr_20.vcf.gz
bcftools view -r 21 final_filtered.vcf.gz -Oz -o by_chr/chr_21.vcf.gz
bcftools view -r 22 final_filtered.vcf.gz -Oz -o by_chr/chr_22.vcf.gz

#Extract genotypes
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