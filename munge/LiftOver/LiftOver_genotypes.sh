
update_build.sh Combo_Barreiro_3projects.BIN.NamesCorrected.171 strand_files/HumanOmni5Exome-4v1-1_A-b38.strand Macrophages_Nedelec_2016

# First generate a frequency file from the chip data
plink -bfile Macrophages_Nedelec_2016.GRCh38 --freq --out Macrophages_Nedelec_2016.GRCh38_plink

# Remove the PLINK format frequency file header
sed -i -e '1d' Macrophages_Nedelec_2016.GRCh38_plink.frq

# Generate a tab-separated frequency file for the major allele
# Include 'chr' tags to the chromosome names in the first 2 columns as this is the notation in hg38
# Add also a SNP column (in format CHR_POS_REF_ALT here columns $1, $4, $6 and $5 of the .bim file) and a header which match to the header format in the panel.frq file
paste <(awk -v OFS='\t' '{print "chr"$1, "chr"$1"_"$4"_"$6"_"$5, $6, $5}' Macrophages_Nedelec_2016.GRCh38.bim) <(awk -v OFS='\t' '{print $5}' Macrophages_Nedelec_2016.GRCh38_plink.frq) | awk 'BEGIN{printf "CHR\tSNP\tREF\tALT\tAF\n"} {print $0}' > Macrophages_Nedelec_2016.GRCh38_major_allele.frq

# Generate a tab-separated frequency file for the minor allele
# Include 'chr' tags to the chromosome names in the first 2 columns as this is the notation in hg38
# Add also a SNP column (in format CHR_POS_REF_ALT here columns $1, $4, $5 and $6 of the .bim file) which match to the format in the panel.frq file
paste <(awk -v OFS='\t' '{print "chr"$1, "chr"$1"_"$4"_"$5"_"$6, $5, $6}' Macrophages_Nedelec_2016.GRCh38.bim) <(awk -v OFS='\t' '{print 1-$5}' Macrophages_Nedelec_2016.GRCh38_plink.frq) > Macrophages_Nedelec_2016.GRCh38_minor_allele.frq

# Concatenate the two files into a single .frq file
cat Macrophages_Nedelec_2016.GRCh38_major_allele.frq Macrophages_Nedelec_2016.GRCh38_minor_allele.frq > Macrophages_Nedelec_2016.GRCh38_vcf_format.frq

# Ensure that chrX uses the GRCh38/hg38 notation
sed -i 's/chr23/chrX/g' Macrophages_Nedelec_2016.GRCh38_vcf_format.frq

#Compare allele frequences against reference
Rscript --no-save ~/software/utils/compare_AF.R Macrophages_Nedelec_2016_vcf_format.frq ~/rocket/datasets/1000G/freq/1000GP_ALL_GRCh38.frq

# Find the alleles that were on - strand and were thus incorrectly flipped with update_build.sh
cat strand_files/HumanOmni5Exome-4v1-1_A-b38.strand | awk '{if ($5 == "-") print $0}' | cut -f 1 > flipped_in_update_build_script.txt

# Flip the alleles back
plink --bfile Macrophages_Nedelec_2016 --flip flipped_in_update_build_script.txt --make-bed --out Macrophages_Nedelec_2016_reflipped

# First generate a frequency file from the chip data
plink -bfile Macrophages_Nedelec_2016_reflipped --freq --out Macrophages_Nedelec_2016_plink


# Check your .bim file for B alleles
grep -P '\tB\t|\tB$' Macrophages_Nedelec_2016.bim | cut -f 2 > Macrophages_Nedelec_2016_Bvariantlist.txt

# Convert the PLINK format to VCF format
plink --bfile Macrophages_Nedelec_2016 --recode vcf-iid --output-chr M --out Macrophages_Nedelec_2016

rs912991 - correct
rs6671356 - flipped


plink --bfile Combo_Barreiro_3projects.BIN.NamesCorrected.171 --recode vcf-iid --output-chr M --out original


#Lift-over
update_build.sh Combo_Barreiro_3projects.BIN.NamesCorrected.171 strand_files/b38/humanomni5exome-4v1_a-b38.strand Macrophages_Nedelec_2016.GRCh38

#Convert to VCF
plink --bfile Macrophages_Nedelec_2016.GRCh38 --recode vcf-iid --output-chr M --out Macrophages_Nedelec_2016.GRCh38

#Add all INFO tags
bcftools +fill-tags Macrophages_Nedelec_2016.GRCh38.vcf -Oz -o Macrophages_Nedelec_2016.GRCh38.AF.vcf.gz

#Filter rare (AC<1) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Macrophages_Nedelec_2016.GRCh38.INFO.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' -Oz -o Macrophages_Nedelec_2016.GRCh38.INFO.filtered.vcf.gz

#Remove variants that do not match to the reference genome
bcftools norm -cx -f ~/rocket/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa Macrophages_Nedelec_2016.GRCh38.INFO.filtered.vcf.gz -Oz -o  Macrophages_Nedelec_2016.GRCh38.INFO.filtered.non_ref.vcf.gz

# Remove duplicates
bcftools norm -d all Macrophages_Nedelec_2016.GRCh38.INFO.filtered.non_ref.vcf.gz -Oz -o test.vcf.gz

# First generate a tab-delimited header for the allele frequency file
echo -e 'CHR\tSNP\tREF\tALT\tAF' > Macrophages_Nedelec_2016.GRCh38.INFO.filtered_chip.frq

# Query the required fields from the VCF file and append to the allele frequency file
bcftools query -f '%CHROM\t%ID\t%REF\t%ALT\t%INFO/AF\n' Macrophages_Nedelec_2016.GRCh38.INFO.filtered.vcf.gz >> Macrophages_Nedelec_2016.GRCh38.INFO.filtered_chip.frq





#GRCh38
#Lift-over
update_build.sh Combo_Barreiro_3projects.BIN.NamesCorrected.171 strand_files/b38/humanomni5exome-4v1_a-b38.strand Macrophages_Nedelec_2016.GRCh38

plink -bfile Macrophages_Nedelec_2016.GRCh38 --reference-allele strand_files/b38/HumanOmni5Exome-4v1_A-b38.strand.RefAlt --make-bed --out Macrophages_Nedelec_2016.GRCh38.RefAlt

plink --bfile Macrophages_Nedelec_2016.GRCh38.RefAlt --recode vcf-iid --output-chr M --out Macrophages_Nedelec_2016.GRCh38.RefAlt

#Fix reference allele
bcftools norm -cs -f ~/rocket/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa Macrophages_Nedelec_2016.GRCh38.RefAlt.vcf -Oz -o  Macrophages_Nedelec_2016.GRCh38.non_ref.vcf.gz

#Remove multi-allelic SNPs
bcftools norm -d all Macrophages_Nedelec_2016.GRCh38.non_ref.vcf.gz | bcftools norm -m+any - | bcftools view -m2 -M2 - -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.vcf.gz 

#Extract test chromosomes
bcftools index Macrophages_Nedelec_2016.GRCh38.biallelic.vcf.gz
bcftools view -r 21 Macrophages_Nedelec_2016.GRCh38.biallelic.vcf.gz -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.vcf.gz

#Rename variants by chr and position

#Add all INFO tags
bcftools +fill-tags Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.vcf.gz | bcftools annotate --set-id +'%CHROM\_%POS' - -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.tags.vcf.gz

#Filter rare (MAF<0.01) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' - -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.filtered.vcf.gz






#GRCh37
#Lift-over
update_build.sh Combo_Barreiro_3projects.BIN.NamesCorrected.171 strand_files/b37/HumanOmni5Exome-4v1_A-b37.strand Macrophages_Nedelec_2016.GRCh38

plink -bfile Macrophages_Nedelec_2016.GRCh38 --reference-allele strand_files/b37/HumanOmni5Exome-4v1_A-b37.strand.RefAlt --make-bed --out Macrophages_Nedelec_2016.GRCh38.RefAlt

plink --bfile Macrophages_Nedelec_2016.GRCh38.RefAlt --recode vcf-iid --output-chr M --out Macrophages_Nedelec_2016.GRCh38.RefAlt

#Fix reference allele
bcftools view -t ^M Macrophages_Nedelec_2016.GRCh38.RefAlt.vcf | bcftools norm -cs -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa - -Oz -o Macrophages_Nedelec_2016.GRCh38.non_ref.vcf.gz

#Remove multi-allelic SNPs
bcftools norm -d all Macrophages_Nedelec_2016.GRCh38.non_ref.vcf.gz | bcftools norm -m+any - | bcftools view -m2 -M2 - -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.vcf.gz 

#Extract test chromosomes
bcftools index Macrophages_Nedelec_2016.GRCh38.biallelic.vcf.gz
bcftools view -r 21 Macrophages_Nedelec_2016.GRCh38.biallelic.vcf.gz -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.vcf.gz

#Rename variants by chr and position

#Add all INFO tags
bcftools +fill-tags Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.vcf.gz | bcftools annotate --set-id +'%CHROM\_%POS' - -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.tags.vcf.gz

#Filter rare (MAF<0.01) and non-HWE varaints and those with abnormal reference alleles
bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' - -Oz -o Macrophages_Nedelec_2016.GRCh38.biallelic.chr21.filtered.vcf.gz


#Recaclulate INFO scores
bcftools +impute-info -Oz -o chr21.imp2.vcf.gz chr21.dose.vcf.gz

#filter based on INFO score
bcftools filter -i 'INFO/INFO > 0.7' chr21.imp2.vcf.gz -Oz -o chr21.filtered.vcf.gz


#Flip non-matched ref and alt alleles
bcftools +fixref chr21_test.vcf.gz -Oz -o chr21_test.flipped.vcf.gz -- -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa -i ~/rocket/datasets/dbSNP/dbSNP_b151_GRCh37p13.vcf.gz

bcftools norm -cs -f ~/rocket/annotations/GRCh37/Homo_sapiens.GRCh37.dna.primary_assembly.fa chr21_test.flipped.vcf.gz | bcftools norm -d all - | bcftools norm -m+any - | bcftools view -m2 -M2 - -Oz -o chr21_test.filtered.vcf.gz


#Check AF distributions
bcftools annotate -c INFO/AF -a ~/rocket/datasets/1000G/GRCh37_allele_frequencies.vcf.gz Macrophages_Nedelec_2016.fixed.vcf.gz | bcftools +af-dist | grep ^PROB > data.dist.txt
