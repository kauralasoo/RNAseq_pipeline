#Convert BED to MAP
plink -bfile ROSMAP_arrayGenotype --recode --tab --out ROSMAP_raw

#Lift over Coordnates to GRCh37
liftOverPlink.py -m ROSMAP_raw.map -p ROSMAP_raw.ped -o ROSMAP_raw_GRCh37 -c NCBI36_to_GRCh37.chain.gz

#Convert back to bed/bim/fam
plink --file ROSMAP_raw_GRCh37 --make-bed --out ROSMAP_raw_GRCh37_bed

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

#Filter final vcf file and add unique variant ids
bcftools filter -i 'MAF[0] > 0.01' Monocytes_Quach_2016_GRCh38.vcf.gz | bcftools annotate --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o Quach_2016_GRCh38.filtered.vcf.gz
bcftools index Quach_2016_GRCh38.filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags Quach_2016_GRCh38.filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\t%R2\n' | gzip > Quach_2016_GRCh38.variant_information.txt.gz
