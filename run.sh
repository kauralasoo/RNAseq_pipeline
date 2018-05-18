
#### GEUVADIS dataset ####
#Run the quantification pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription.snakefile processed/GEUVADIS/out.txt --jobs 100 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete

snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/GEUVADIS/out.txt --jobs 100 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete


#### TwinsUK dataset ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE.snakefile processed/TwinsUK/out.txt --jobs 1 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

#Run the alignment pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/TwinsUK/out.txt --jobs 100 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/TwinsUK/out.txt --jobs 100 --configfile configs/TwinsUK_config.yaml --rerun-incomplete



### Nedelec 2017 dataset ####
#Convert .sra files to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s sra_to_fastq.snakefile processed/Macrophages_Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_config.yaml --rerun-incomplete

#Quantify transcription
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_SE.snakefile processed/Macrophages_Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_config.yaml --rerun-incomplete
