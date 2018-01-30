
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
