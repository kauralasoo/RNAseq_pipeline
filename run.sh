
#### GEUVADIS dataset ####
#Run the quantification pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription.snakefile processed/GEUVADIS/out.txt --jobs 100 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete

snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/GEUVADIS/out.txt --jobs 20 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete


#### TwinsUK dataset ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE.snakefile processed/TwinsUK/out.txt --jobs 1 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

#Run the alignment pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/TwinsUK/out.txt --jobs 100 --configfile configs/TwinsUK_config.yaml --rerun-incomplete


#### Nedelec 2016 dataset ####
#Convert .sra files to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s sra_to_fastq.snakefile processed/Macrophages_Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_config.yaml --rerun-incomplete

#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Macrophages_Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_config.yaml --rerun-incomplete


#### Quach 2016 dataset ####
#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Monocytes_Quach_2016/out.txt --jobs 1 --configfile configs/Quach_2016_config.yaml --rerun-incomplete


#### BLUEPRINT dataset ####
#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/BLUEPRINT/out.txt --jobs 1 --configfile configs/BLUEPRINT_SE_config.yaml --rerun-incomplete
snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription_PE.snakefile processed/BLUEPRINT/out.txt --jobs 20 --configfile configs/BLUEPRINT_PE_config.yaml --rerun-incomplete

snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription_PE.snakefile processed/BLUEPRINT/out.txt --jobs 20 --configfile configs/BLUEPRINT_all_config.yaml --rerun-incomplete


#### Alasoo et al ####
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/Macrophages_Alasoo_2018/out.txt --jobs 20 --configfile configs/Alasoo_2018_config.yaml --rerun-incomplete

#### Fairfax et al ####
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Fairfax/out.txt --jobs 1 --configfile configs/Fairfax_config.yaml --rerun-incomplete