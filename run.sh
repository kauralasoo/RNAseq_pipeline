#### GEUVADIS dataset ####
#Run the quantification pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription.snakefile processed/GEUVADIS/out.txt --jobs 100 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete

snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/GEUVADIS/out.txt --jobs 20 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete

#Map QTLs
snakemake --cluster scripts/snakemake_submit_UT.py -np -s map_QTLs.snakefile processed/GEUVADIS/out.txt --configfile configs/GEUVADIS_config.yaml --rerun-incomplete --jobs 100


#### TwinsUK dataset ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE.snakefile processed/TwinsUK/out.txt --jobs 1 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

#Run the alignment pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/TwinsUK/out.txt --jobs 40 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes_no_R2.snakefile -np processed/TwinsUK/genotypes/GRCh38/TwinsUK_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

#Map QTLs
snakemake --cluster scripts/snakemake_submit_UT.py -np -s map_QTLs.snakefile processed/TwinsUK/out.txt --configfile configs/TwinsUK_config.yaml --rerun-incomplete --jobs 100


#### Nedelec 2016 dataset ####
#Convert .sra files to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s sra_to_fastq.snakefile processed/Macrophages_Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_config.yaml --rerun-incomplete

#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Macrophages_Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_config.yaml --rerun-incomplete

#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -p processed/Macrophages_Nedelec_2016/genotypes/GRCh38/Macrophages_Nedelec_2016_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 10 --rerun-incomplete


#### Quach 2016 dataset ####
#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Monocytes_Quach_2016/out.txt --jobs 20 --configfile configs/Quach_2016_config.yaml --rerun-incomplete

#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -p processed/Monocytes_Quach_2016/genotypes/GRCh38/Monocytes_Quach_2016_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 10 --rerun-incomplete

#Map QTLs
snakemake --cluster scripts/snakemake_submit_UT.py -np -s map_QTLs.snakefile processed/Monocytes_Quach_2016/out.txt --configfile configs/Quach_2016_config.yaml --rerun-incomplete --jobs 200


#### BLUEPRINT dataset ####
#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/BLUEPRINT/out.txt --jobs 30 --configfile configs/BLUEPRINT_SE_config.yaml --rerun-incomplete
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/BLUEPRINT/out.txt --jobs 35 --configfile configs/BLUEPRINT_PE_config.yaml --rerun-incomplete

#Run using all sample names at the same time
snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription_PE.snakefile processed/BLUEPRINT/out.txt --jobs 20 --configfile configs/BLUEPRINT_all_config.yaml --rerun-incomplete


#### Alasoo et al ####
snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription_PE.snakefile processed/Macrophages_Alasoo_2018/out.txt --jobs 8 --configfile configs/Alasoo_2018_config.yaml --rerun-incomplete

#### Fairfax et al ####
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Fairfax/out.txt --jobs 1 --configfile configs/Fairfax_config.yaml --rerun-incomplete




#### GENCORD ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE_tryHPC.snakefile processed/GENCORD/out.txt --jobs 20 --configfile configs/GENCORD_config.yaml --rerun-incomplete

#Run the alignment pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/GENCORD/out.txt --jobs 30 --configfile configs/GENCORD_config_align.yaml --rerun-incomplete

#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/GENCORD/out.txt --configfile configs/CrossMap_config.yaml --jobs 10 --rerun-incomplete

#Map QTLs
snakemake --cluster scripts/snakemake_submit_UT.py -np -s map_QTLs.snakefile processed/GENCORD/out.txt --configfile configs/GENCORD_config_align.yaml --rerun-incomplete --jobs 100



#### Schwartzentruber ####
#Convert crams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s cram_to_fastq_PE.snakefile processed/SensoryNeurons/out.txt --jobs 20 --configfile configs/Schwartzentruber_2018_config_cram_to_fastq.yaml --rerun-incomplete

#Run the alignment pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/Schwartzentruber_2018/out.txt --jobs 30 --configfile configs/Schwartzentruber_2018_config_align.yaml --rerun-incomplete



#### HipSci ####
#Convert crams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s cram_to_fastq_PE.snakefile processed/HipSci/out.txt --jobs 30 --configfile configs/HipSci_config_cram_to_fastq.yaml --rerun-incomplete

#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/HipSci/genotypes/GRCh38/HipSci_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs --rerun-incomplete



##### CEDAR ####
#CrossMap genotypes
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -np processed/CEDAR/out.txt --configfile configs/CrossMap_config.yaml --jobs 1 --rerun-incomplete




#Run Samtools indexing on HPC
snakemake --cluster scripts/snakemake_submit_UT.py -np -s samtools_index.snakefile processed/BLUEPRINT/out.txt --jobs 20 --configfile configs/BLUEPRINT_all_config.yaml --rerun-incomplete

#Run Samtools MBV analysis on HPC
snakemake --cluster scripts/snakemake_submit_UT.py -np -s mbv_analysis.snakefile processed/BLUEPRINT/out.txt --jobs 20 --configfile configs/BLUEPRINT_all_config.yaml --rerun-incomplete


#### van_de_Bunt_2015 ####
snakemake --cluster scripts/snakemake_submit_UT.py -s CrossMap_genotypes.snakefile -p processed/van_de_Bunt_2015/genotypes/GRCh38/van_de_Bunt_2015_GRCh38.vcf.gz --configfile configs/CrossMap_config.yaml --jobs 22 --rerun-incomplete

