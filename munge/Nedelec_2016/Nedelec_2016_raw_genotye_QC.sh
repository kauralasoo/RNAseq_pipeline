module load bcftools-1.8

#Update build
module load plink-1.9.0
update_build.sh Combo_Barreiro_3projects.BIN.NamesCorrected.171 /gpfs/hpchome/a72094/rocket/datasets/imputation/strand_files/b37/HumanOmni5Exome-4v1_A-b37.strand Nedelec_2016

#Fix Ref and Alt alleles with PLINK
plink -bfile Nedelec_2016 --reference-allele /gpfs/hpchome/a72094/rocket/datasets/imputation/strand_files/b37/HumanOmni5Exome-4v1_A-b37.strand.RefAlt --make-bed --out Nedelec_2016_RefAlt

#Convert to VCF
plink --bfile Nedelec_2016_RefAlt --recode vcf-iid --out Nedelec_2016_RefAlt
