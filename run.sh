snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/GEUVADIS/STAR/ERR188374/ERR188374.Aligned.sortedByCoord.out.bam --jobs 1 --configfile configs/GEUVADIS_config.yaml

snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/GEUVADIS/bigwig/ERR188307.str1.bw --jobs 1 --configfile configs/GEUVADIS_config.yaml


"processed/{study}/salmon/{annotation}/{sample}/quant.sf"


snakemake --cluster scripts/snakemake_submit_UT.py -p -s quantify_transcription.snakefile processed/GEUVADIS/out.txt --jobs 100 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete

snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/GEUVADIS/out.txt --jobs 100 --configfile configs/GEUVADIS_config.yaml --rerun-incomplete


rm /tmp/*.fifo*


#### TwinsUK dataset ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE.snakefile processed/TwinsUK/out.txt --jobs 1 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

#Run the alignment pipeline
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription.snakefile processed/TwinsUK/out.txt --jobs 100 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

