#### ImmVar dataset ####
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/ImmVar/genotypes/GRCh38/ImmVar_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

#### Schmiedel_2018 dataset ####
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/Schmiedel_2018/genotypes/GRCh38/Schmiedel_2018_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

#### Kasela_2017 dataset ####
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/Kasela_2017/genotypes/GRCh38/Kasela_2017_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

#### van_de_Bunt_2015 ####
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -p processed/van_de_Bunt_2015/genotypes/GRCh38/van_de_Bunt_2015_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

##### Fairfax_2014 ####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/Fairfax_2014/genotypes/GRCh38/Fairfax_2014_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 1 --rerun-incomplete

#### Garieri 2017 ####
#Extract samples from 1KG VCFs
snakemake -p --snakefile extract_genotypes.snakefile --cluster ../../scripts/snakemake_submit_UT.py --jobs 10

##### CEDAR ####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/CEDAR/out.txt --configfile configs/CrossMap_config.yaml --jobs 1 --rerun-incomplete

##### HipSci #####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/HipSci/genotypes/GRCh38/HipSci_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs --rerun-incomplete

##### GENCORD #####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/GENCORD/out.txt --configfile configs/CrossMap_config.yaml --jobs 10 --rerun-incomplete

##### TwinsUK #####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes_no_R2.snakefile -np processed/TwinsUK/genotypes/GRCh38/TwinsUK_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

##### Nedelec_2016 #####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -p processed/Nedelec_2016/genotypes/GRCh38/Macrophages_Nedelec_2016_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 10 --rerun-incomplete

#### Quach 2016 dataset ####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -p processed/Monocytes_Quach_2016/genotypes/GRCh38/Monocytes_Quach_2016_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 10 --rerun-incomplete
