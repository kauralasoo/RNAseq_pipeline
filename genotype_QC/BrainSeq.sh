#Update the .fam file
cp ~/datasets/controlled_access/SampleArcheology/studies/BrainSeq/LIBD_Brain_DLPFC_szControls_imputed.fam .

#Convert PLINK to VCF
plink --bfile LIBD_Brain_DLPFC_szControls_imputed --recode vcf-iid --out BrainSeq_full

##### HAP650 ####
wget https://www.well.ox.ac.uk/~wrayner/strand/HumanHap650Yv3_A-b37-strand.zip

#Extract samples
bcftools view -S ~/datasets/controlled_access/SampleArcheology/studies/BrainSeq/h650_chip_samples.txt BrainSeq_full.vcf.gz -Oz -o BrainSeq_h650.vcf.gz

#Extract h650 positions
cut -f 2,3 HumanHap650Yv3_A-b37.strand  > h650_positions.txt
bcftools view -R h650_positions.txt BrainSeq_h650.vcf.gz -Oz -o BrainSeq_h650_filtered.vcf.gz

#Convert back to plink format for QC
plink --vcf BrainSeq_h650_filtered.vcf.gz --make-bed --out BrainSeq_h650


##### 1M ####
wget https://www.well.ox.ac.uk/~wrayner/strand/Human1MDuov3_C-b37-strand.zip

#Extract samples
bcftools view -S /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/studies/BrainSeq/1M_chip_samples.txt BrainSeq_full.vcf.gz -Oz -o BrainSeq_1M.vcf.gz

#Extract h650 positions
cut -f 2,3 Human1M-Duov3_C-b37.strand  > 1m_positions.txt
bcftools view -R 1m_positions.txt BrainSeq_1M.vcf.gz -Oz -o BrainSeq_1M_filtered.vcf.gz

#Convert back to plink format for QC
plink --vcf BrainSeq_1M_filtered.vcf.gz --make-bed --out BrainSeq_1M


#Merge imputation results
bcftools index ../1M/GRCh38/BrainSeq_1M.vcf.gz 
bcftools index ../h650/GRCh38/BrainSeq_h650.vcf.gz
bcftools merge ../1M/GRCh38/BrainSeq_1M.vcf.gz ../h650/GRCh38/BrainSeq_h650.vcf.gz -Oz -o BrainSeq_GRCh38_merged.vcf.gz

#Filter missing
bcftools filter -i 'F_MISSING < 0.05 & MAF[0] > 0.01' BrainSeq_GRCh38_merged.vcf.gz -Oz -o BrainSeq_GRCh38_filtered.vcf.gz
