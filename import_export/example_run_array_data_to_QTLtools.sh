Rscript array_data_to_QTLtools_input.R\
 -g "../metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz"\
 -s "../metadata/cleaned/Fairfax_2014.tsv"\
 -e "/gpfs/hpc/home/a72094/projects/RNAseq_pipeline/results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz"\
 -v "/gpfs/hpchome/a72094/hpc/datasets/controlled_access/Fairfax_2014/genotypes/Michigan_GRCh37_1KGPhase3_061118/GRCh38/Fairfax_2014_GRCh38.variant_information.txt.gz"\
 --qtlutils "../../eQTLUtils"\
 -o "../processed/Fairfax_2014/qtltools/input/array/" #-c 1000001 -m 6