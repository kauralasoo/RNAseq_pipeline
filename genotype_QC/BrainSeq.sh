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
