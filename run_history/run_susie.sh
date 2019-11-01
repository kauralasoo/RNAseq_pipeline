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