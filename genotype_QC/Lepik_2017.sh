#Split multi-allelic varaints
bcftools norm -m -any EGCUT_MAF001.vcf.gz -Oz -o EGCUT_MAF001_no_multiallelic.vcf.gz

#Extract chromosomes from the VCF file
bcftools view -r 1 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr1.dose.vcf.gz
bcftools view -r 2 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr2.dose.vcf.gz
bcftools view -r 3 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr3.dose.vcf.gz
bcftools view -r 4 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr4.dose.vcf.gz
bcftools view -r 5 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr5.dose.vcf.gz
bcftools view -r 6 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr6.dose.vcf.gz
bcftools view -r 7 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr7.dose.vcf.gz
bcftools view -r 8 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr8.dose.vcf.gz
bcftools view -r 9 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr9.dose.vcf.gz
bcftools view -r 10 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr10.dose.vcf.gz
bcftools view -r 11 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr11.dose.vcf.gz
bcftools view -r 12 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr12.dose.vcf.gz
bcftools view -r 13 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr13.dose.vcf.gz
bcftools view -r 14 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr14.dose.vcf.gz
bcftools view -r 15 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr15.dose.vcf.gz
bcftools view -r 16 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr16.dose.vcf.gz
bcftools view -r 17 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr17.dose.vcf.gz
bcftools view -r 19 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr19.dose.vcf.gz
bcftools view -r 20 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr20.dose.vcf.gz
bcftools view -r 21 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr21.dose.vcf.gz
bcftools view -r 22 EGCUT_MAF001_no_multiallelic.vcf.gz -Oz -o GRCh37/chr22.dose.vcf.gz


