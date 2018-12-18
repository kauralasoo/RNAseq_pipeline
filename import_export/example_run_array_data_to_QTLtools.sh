Rscript array_data_to_QTLtools_in.R\
 -g "../metadata/gene_metadata/HumanHT-12_V4_gene_metadata.txt.gz"\
 -s "../metadata/cleaned/Fairfax_2014.tsv"\
 -e "../results/expression_matrices/HumanHT-12_V4/Fairfax_2014.tsv.gz"\
 -v "../../temp/Fairfax_2014_GRCh38.variant_information.txt.gz"\
 --qtlutils "../../eQTLUtils"\
 -o "../processed/Fairfax_2014/qtltools/input/array/" #-c 1000001 -m 6