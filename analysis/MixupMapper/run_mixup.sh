#Convert genotypes to TRITYPER format
java -jar ~/software/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i CEDAR_GRCh38.filtered.renamed -I VCF -o CEDAR -O TRITYPER

#Run MixupMapper
java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out mixup_out\
 --inexp monocytes_experession.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --snps snp_list.txt\
 --gte geno_pheno_coupling.tsv
