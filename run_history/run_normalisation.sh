### CEDAR ###
./nextflow run normalisation.nf -profile tartu_hpc \
    -resume \
    --study_name CEDAR\
    --is_microarray\
    --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/CEDAR.tsv.gz\
    --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/CEDAR.tsv\
    --vcf_file /gpfs/hpc/projects/genomic_references/CEDAR/genotypes/Michigan_GRCh37_1KGPhase3_220918/GRCh38/CEDAR_GRCh38.filtered.renamed.vcf.gz\
    --outdir HumanHT-12_V4

### Fairfax_2012 ###
./nextflow run normalisation.nf -profile tartu_hpc \
    -resume \
    --study_name Fairfax_2012\
    --is_microarray\
    --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/Fairfax_2012.tsv.gz\
    --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Fairfax_2012.tsv\
    --vcf_file /gpfs/hpc/projects/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
    --outdir HumanHT-12_V4

### Fairfax_2014 ###
./nextflow run normalisation.nf -profile tartu_hpc \
    -resume \
    --study_name Fairfax_2014\
    --is_microarray\
    --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/Fairfax_2014.tsv.gz\
    --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Fairfax_2014.tsv\
    --vcf_file /gpfs/hpc/projects/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
    --outdir HumanHT-12_V4

### Naranbhai_2015 ###
./nextflow run normalisation.nf -profile tartu_hpc \
    -resume \
    --study_name Naranbhai_2015\
    --is_microarray\
    --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/Naranbhai_2015.tsv.gz\
    --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Naranbhai_2015.tsv\
    --vcf_file /gpfs/hpc/projects/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.filtered.renamed.vcf.gz\
    --outdir HumanHT-12_V4

### Kasela_2017 ###
./nextflow run normalisation.nf -profile tartu_hpc \
    -resume \
    --study_name Kasela_2017\
    --is_microarray\
    --microarray_exp_matrix_path /gpfs/hpc/projects/eQTLCatalogue/processed/expression_matrices/HumanHT-12_V4/raw/Kasela_2017.tsv.gz\
    --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Kasela_2017.tsv\
    --vcf_file /gpfs/hpc/projects/EGCUT_eQTLs/Kasela_2017/genotypes/Michigan_GRCh37_1KGPhase3_220119/GRCh38/Kasela_2017_GRCh38.filtered.vcf.gz\
    --outdir HumanHT-12_V4



##### Normalise RNA-seq read counts for Kateryna

#GENCORD
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name GENCORD\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/GENCORD/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/GENCORD.tsv\
 --vcf_file /gpfs/hpc/projects/GENCORD/genotypes/Michigan_GRCh37_1KGPhase3_220918/GRCh38/GENCORD_GRCh38.filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#Alasoo_2018
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Alasoo_2018\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/Alasoo_2018/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Alasoo_2018.tsv\
 --vcf_file /gpfs/hpc/projects/HipSci/Alasoo_2018/genotypes/Alasoo_2018_GRCh38.filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#Lepik_2017
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Lepik_2017\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/Lepik_2017/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Lepik_2017.tsv\
 --vcf_file /gpfs/hpc/projects/EGCUT_eQTLs/Lepik_2017/genotypes/GRCh38/Lepik_2017_GRCh38.filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#HipSci
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name HipSci\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/HipSci/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/HipSci.tsv\
 --vcf_file /gpfs/hpc/projects/HipSci/genotypes/HipSci_GRCh38.filtered.needed.323.samples.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#Nedelec_2016
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Nedelec_2016\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/Nedelec_2016/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Nedelec_2016.tsv\
 --vcf_file /gpfs/hpc/projects/Nedelec_2016/genotypes/Michigan_GRCh37_130618/GRCh38/Nedelec_2016_GRCh38.filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#Quach_2016
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Quach_2016\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/Quach_2016/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Quach_2016.tsv\
 --vcf_file /gpfs/hpc/projects/Quach_2016/genotypes/Michigan_GRCh37_Phase3_240918/GRCh38/Quach_2016_GRCh38.filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

 #GEUVADIS
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name GEUVADIS\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/GEUVADIS/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/GEUVADIS.tsv\
 --vcf_file /gpfs/hpc/projects/genomic_references/GEUVADIS/genotypes/GEUVADIS_GRCh38_filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#FUSION
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name FUSION\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/FUSION/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/FUSION.tsv\
 --vcf_file /gpfs/hpc/projects/FUSION/genotypes/FUSION_illumina.MAF001.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

 #Schmiedel_2018
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Schmiedel_2018\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/Schmiedel_2018/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Schmiedel_2018.tsv\
 --vcf_file /gpfs/hpc/projects/Schmiedel_2018/genotypes/Michigan_GRCh37_Phase3_210819/GRCh38/Schmiedel_2018.MAF001.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#ROSMAP
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name ROSMAP\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/ROSMAP/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/ROSMAP.tsv\
 --vcf_file /gpfs/hpc/projects/ROSMAP/genotypes/Michigan_GRCh37_Phase3_200819/merged/ROSMAP_GRCh38_filtered.no_multi.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

 #BrainSeq
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name BrainSeq\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/BrainSeq/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/BrainSeq.tsv\
 --vcf_file /gpfs/hpc/projects/BrainSeq/genotypes/Michigan_GRCh37_Phase3_250819/merged/BrainSeq_GRCh38_filtered.no_multi.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#Schwartzentruber_2018
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name Schwartzentruber_2018\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/Schwartzentruber_2018/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/Schwartzentruber_2018.tsv\
 --vcf_file /gpfs/hpc/projects/HipSci/Schwartzentruber_2018/Schwartzentruber_2018_needed_samples.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

#van_de_Bunt_2015
./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name van_de_Bunt_2015\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/van_de_Bunt_2015/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/van_de_Bunt_2015.tsv\
 --vcf_file /gpfs/hpc/projects/van_de_Bunt_2015/genotypes/Michigan_GRCh37_1KGPhase3_051018/GRCh38/van_de_Bunt_2015_GRCh38.filtered.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd


 #BLUEPRINT_SE
 ./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name BLUEPRINT_SE\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/BLUEPRINT_SE/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/BLUEPRINT_SE.tsv\
 --vcf_file  /gpfs/hpc/home/a72094/hpcproject/BLUEPRINT/genotypes/BLUEPRINT_v2/GRCh38/BLUEPRINT.MAF001.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

 #BLUEPRINT_PE
 ./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name BLUEPRINT_PE\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/BLUEPRINT_PE/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/BLUEPRINT_PE.tsv\
 --vcf_file  /gpfs/hpc/home/a72094/hpcproject/BLUEPRINT/genotypes/BLUEPRINT_v2/GRCh38/BLUEPRINT.MAF001.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd

  #BLUEPRINT_PE
 ./nextflow run normalisation.nf -profile tartu_hpc -resume\
 --study_name TwinsUK\
 --quant_results_path /gpfs/hpc/projects/eQTLCatalogue/rnaseq/TwinsUK/\
 --sample_meta_path /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/studies/cleaned/TwinsUK.tsv\
 --vcf_file  /gpfs/hpc/projects/TwinsUK/genotypes/vcf/1KG_Phase3_200520/merged/TwinsUK.MAF001.vcf.gz\
 --skip_exon_norm\
 --skip_tx_norm\
 --skip_leafcutter_norm\
 --skip_txrev_norm\
 --outdir RNAseq\
 -process.queue amd