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