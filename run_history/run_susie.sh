# HumanHT12V4
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/HumanHT12V4.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

 # Kolberg_2020
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/Kolberg_2020.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

 # Fairfax_2014
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/finemapping/Fairfax_2014.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

 # Fairfax_2018
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/Fairfax_2018.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

#Raj_2014
NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/Raj_2014.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

  # Schmiedel_2018
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/RNAseq_DS.tsv\
 --cisdistance 1000000\
 --n_batches 400\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

# BLUEPRINT
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/RNAseq_BLUEPRINT.tsv\
 --cisdistance 1000000\
 --n_batches 400\
 --permuted 'false'\
 --vcf_genotype_field GT\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

# Other dosage datasets
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/RNAseq_DS_other.tsv\
 --cisdistance 1000000\
 --n_batches 400\
 --permuted 'false'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200