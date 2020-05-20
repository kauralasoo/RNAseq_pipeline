#Concat all chromosomes
bcftools concat chr1.dose.vcf.gz chr2.dose.vcf.gz chr3.dose.vcf.gz chr4.dose.vcf.gz chr5.dose.vcf.gz chr6.dose.vcf.gz chr7.dose.vcf.gz chr8.dose.vcf.gz chr9.dose.vcf.gz chr10.dose.vcf.gz chr11.dose.vcf.gz chr12.dose.vcf.gz chr13.dose.vcf.gz chr14.dose.vcf.gz chr15.dose.vcf.gz chr16.dose.vcf.gz chr17.dose.vcf.gz chr18.dose.vcf.gz chr19.dose.vcf.gz chr20.dose.vcf.gz chr21.dose.vcf.gz chr22.dose.vcf.gz -Oz -o TwinsUK_all.vcf.gz

#Extract samples for the two arrays
bcftools view -S ~/eQTLCatalogue/SampleArcheology/studies/TwinsUK/array1_samples.txt TwinsUK_all.vcf.gz -Oz -o TwinsUK_Human610.vcf.gz
bcftools view -S ~/eQTLCatalogue/SampleArcheology/studies/TwinsUK/array2_samples.txt TwinsUK_all.vcf.gz -Oz -o TwinsUK_HumanHap300.vcf.gz

#Exptract SNP positions from the two strand files
cut -f2,3 strand_files/HumanHap300v2_A-b37.strand > HumanHap300_positions.tsv
cut -f2,3 strand_files/Human610-Quadv1_B-b37.strand > Human610_positions.tsv

#Keep the correct SNPs only
bcftools index TwinsUK_Human610.vcf.gz
bcftools view -R Human610_positions.tsv TwinsUK_Human610.vcf.gz -Oz -o TwinsUK_Human610.array_variants.vcf.gz

bcftools index TwinsUK_HumanHap300.vcf.gz
bcftools view -R HumanHap300_positions.tsv TwinsUK_HumanHap300.vcf.gz -Oz -o TwinsUK_HumanHap300.array_variants.vcf.gz

#Fill in imputation quality score
bcftools +impute-info TwinsUK_Human610.array_variants.vcf.gz -Oz -o TwinsUK_Human610.array_variants.INFO.vcf.gz
bcftools +impute-info TwinsUK_HumanHap300.array_variants.vcf.gz -Oz -o TwinsUK_HumanHap300.array_variants.INFO.vcf.gz

#Keep variant with very high quality scores
bcftools filter -i 'INFO/INFO = 1' TwinsUK_Human610.array_variants.INFO.vcf.gz -Oz -o TwinsUK_Human610.array_variants.INFO_filtered.vcf.gz
bcftools filter -i 'INFO/INFO = 1' TwinsUK_HumanHap300.array_variants.INFO.vcf.gz -Oz -o TwinsUK_HumanHap300.array_variants.INFO_filtered.vcf.gz

#Sort by position
bcftools sort TwinsUK_HumanHap300.array_variants.INFO_filtered.vcf.gz -Oz -o TwinsUK_HumanHap300.genotyped.vcf.gz
bcftools sort TwinsUK_Human610.array_variants.INFO_filtered.vcf.gz -Oz -o TwinsUK_Human610.genotyped.vcf.gz

#Split by chromosome
bcftools view -r 1 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_1.vcf.gz
bcftools view -r 2 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_2.vcf.gz
bcftools view -r 3 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_3.vcf.gz
bcftools view -r 4 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_4.vcf.gz
bcftools view -r 5 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_5.vcf.gz
bcftools view -r 6 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_6.vcf.gz
bcftools view -r 7 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_7.vcf.gz
bcftools view -r 8 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_8.vcf.gz
bcftools view -r 9 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_9.vcf.gz
bcftools view -r 10 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_10.vcf.gz
bcftools view -r 11 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_11.vcf.gz
bcftools view -r 12 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_12.vcf.gz
bcftools view -r 13 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_13.vcf.gz
bcftools view -r 14 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_14.vcf.gz
bcftools view -r 15 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_15.vcf.gz
bcftools view -r 16 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_16.vcf.gz
bcftools view -r 17 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_17.vcf.gz
bcftools view -r 18 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_18.vcf.gz
bcftools view -r 19 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_19.vcf.gz
bcftools view -r 20 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_20.vcf.gz
bcftools view -r 21 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_21.vcf.gz
bcftools view -r 22 TwinsUK_HumanHap300.genotyped.vcf.gz -Oz -o by_chr/chr_22.vcf.gz

bcftools view -r 1 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_1.vcf.gz
bcftools view -r 2 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_2.vcf.gz
bcftools view -r 3 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_3.vcf.gz
bcftools view -r 4 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_4.vcf.gz
bcftools view -r 5 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_5.vcf.gz
bcftools view -r 6 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_6.vcf.gz
bcftools view -r 7 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_7.vcf.gz
bcftools view -r 8 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_8.vcf.gz
bcftools view -r 9 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_9.vcf.gz
bcftools view -r 10 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_10.vcf.gz
bcftools view -r 11 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_11.vcf.gz
bcftools view -r 12 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_12.vcf.gz
bcftools view -r 13 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_13.vcf.gz
bcftools view -r 14 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_14.vcf.gz
bcftools view -r 15 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_15.vcf.gz
bcftools view -r 16 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_16.vcf.gz
bcftools view -r 17 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_17.vcf.gz
bcftools view -r 18 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_18.vcf.gz
bcftools view -r 19 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_19.vcf.gz
bcftools view -r 20 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_20.vcf.gz
bcftools view -r 21 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_21.vcf.gz
bcftools view -r 22 TwinsUK_Human610.genotyped.vcf.gz -Oz -o by_chr/chr_22.vcf.gz