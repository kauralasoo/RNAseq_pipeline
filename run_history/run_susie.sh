 ##### Kolberg_2020 ####
nextflow run main.nf \
 --qtl_results /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/finemapping/Kolberg_2020.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'true'\
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
 -resume\
 -profile finemapping\
 -executor.queueSize 1

# HumanHT12V4
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/finemapping/HumanHT12V4.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
 --vcf_genotype_field DS\
 -resume\
 -profile finemapping\
 -executor.queueSize 200

 # RNA-seq gene counts
 NXF_VER=18.10.1 nextflow run main.nf \
 --qtl_results /gpfs/hpc/home/a72094/datasets/controlled_access/SampleArcheology/finemapping/RNAseq_gene_counts.tsv\
 --cisdistance 1000000\
 --n_batches 200\
 --permuted 'false'\
 --eqtlutils '/gpfs/hpc/home/a72094/projects/eQTLUtils'\
 -resume\
 -profile finemapping\
 -executor.queueSize 200