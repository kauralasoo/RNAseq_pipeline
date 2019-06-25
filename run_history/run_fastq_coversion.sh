#### TwinsUK dataset ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE.snakefile processed/TwinsUK/out.txt --jobs 1 --configfile configs/TwinsUK_config.yaml --rerun-incomplete

#### Nedelec 2016 dataset ####
#Convert .sra files to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s sra_to_fastq.snakefile processed/Nedelec_2016/out.txt --jobs 1 --configfile configs/Nedelec_2016_config.yaml --rerun-incomplete

#### GENCORD ####
#Convert bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE_tryHPC.snakefile processed/GENCORD/out.txt --jobs 20 --configfile configs/GENCORD_config.yaml --rerun-incomplete

#### Schwartzentruber ####
#Convert crams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s cram_to_fastq_PE.snakefile processed/SensoryNeurons/out.txt --jobs 20 --configfile configs/Schwartzentruber_2018_config_cram_to_fastq.yaml --rerun-incomplete

#### HipSci ####
#Convert crams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s cram_to_fastq_PE.snakefile processed/HipSci/out.txt --jobs 30 --configfile configs/HipSci_config_cram_to_fastq.yaml --rerun-incomplete

#### van_de_Bunt_2015 ####
# Convert Bams to fastq
snakemake --cluster scripts/snakemake_submit_UT.py -np -s bam_to_fastq_PE_filename.snakefile processed/Bunt_2015/out.txt --jobs 30 --configfile configs/Bunt_2015_config_bam_to_fastq.yaml --rerun-incomplete

#### Ye_2018 dataset ####
#Convert .sra files to fastq
snakemake --cluster ~/hpc/projects/RNAseq_pipeline/scripts/snakemake_submit_UT.py -p -s dbGaP_sra_to_fastq.snakefile out.txt --jobs 50 --configfile ~/hpc/projects/RNAseq_pipeline/configs/Ye_2018_config.yaml --rerun-incomplete

#### Lepik_2017 dataset ####
snakemake --cluster scripts/snakemake_submit_UT.py -np -s quantify_transcription_PE.snakefile processed/Lepik_2017/out.txt --jobs 30 --configfile configs/Lepik_2017_config.yaml --rerun-incomplete

### Run coloc all QTLs agains all GWAS traits
snakemake --cluster scripts/snakemake_submit_UT.py -np -s run_coloc.snakefile results/coloc/coloc_out.txt --jobs 300 --configfile configs/coloc_config.yaml --rerun-incomplete

### Extract variant effects
Rscript scripts/extract_variants_from_summaries.R -p results/extracted_variants/lead_vars.txt -q results/extracted_variants/qtl_groups.txt -o results/extracted_variants/results.txt
Rscript scripts/extract_variants_from_summaries.R -p metadata/coloc/featureCounts_coloc_examples.txt -q results/extracted_variants/qtl_groups.txt -o results/extracted_variants/featureCounts_results.txt
