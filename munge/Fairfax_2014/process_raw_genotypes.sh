#Extract sample ids from genotypes
zcat Oxford_eQTL_cohort_genotype2.txt.gz | tail -n+2 | cut -f 2 | uniq > sample_names_2.txt

