module load bcftools-1.8

#Update build
module load plink-1.9.0
update_build.sh Combo_Barreiro_3projects.BIN.NamesCorrected.171 /gpfs/hpchome/a72094/rocket/datasets/imputation/strand_files/b37/HumanOmni5Exome-4v1_A-b37.strand Nedelec_2016

#Fix Ref and Alt alleles with PLINK
plink -bfile Nedelec_2016 --reference-allele /gpfs/hpchome/a72094/rocket/datasets/imputation/strand_files/b37/HumanOmni5Exome-4v1_A-b37.strand.RefAlt --make-bed --out Nedelec_2016_RefAlt

#Convert to VCF
plink --bfile Nedelec_2016_RefAlt --recode vcf-iid --out Nedelec_2016_RefAlt

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.01' Macrophages_Nedelec_2016_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Nedelec_2016_GRCh38.filtered.vcf.gz
bcftools index Nedelec_2016_GRCh38.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags Nedelec_2016_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\t%R2\n' | gzip > Nedelec_2016_GRCh38.variant_information.txt.gz
