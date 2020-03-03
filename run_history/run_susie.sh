# HumanHT12V4
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/projects/eQTLCatalogue/SampleArcheology/finemapping/HumanHT12V4.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
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
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
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

  # RNA-seq dosage field
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