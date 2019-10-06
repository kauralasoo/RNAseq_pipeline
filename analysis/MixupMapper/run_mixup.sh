#Convert genotypes to TRITYPER format
java -jar ~/software/GenotypeHarmonizer-1.4.20-SNAPSHOT/GenotypeHarmonizer.jar -i CEDAR_GRCh38.filtered.renamed -I VCF -o CEDAR -O TRITYPER

#Run MixupMapper
java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out monocyte_CD14\
 --inexp monocyte_CD14.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte monocyte_CD14.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

 java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out B-cell_CD19\
 --inexp B-cell_CD19.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte B-cell_CD19.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

  java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out T-cell_CD8\
 --inexp T-cell_CD8.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte T-cell_CD8.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

  java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out T-cell_CD4\
 --inexp T-cell_CD4.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte T-cell_CD4.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

  java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out neutrophil_CD15\
 --inexp neutrophil_CD15.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte neutrophil_CD15.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

  java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out platelet\
 --inexp platelet.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte platelet.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

  java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out transverse_colon\
 --inexp transverse_colon.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte transverse_colon.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

  java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out ileum\
 --inexp ileum.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte ileum.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

   java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in CEDAR\
 --out rectum\
 --inexp rectum.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte rectum.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

#Fairfax_2014
java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in Fairfax_2014\
 --out monocyte_naive\
 --inexp monocyte_naive.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte monocyte_naive.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

 java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in Fairfax_2014\
 --out monocyte_LPS24\
 --inexp monocyte_LPS24.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte monocyte_LPS24.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in Fairfax_2014\
 --out monocyte_IFN24\
 --inexp monocyte_IFN24.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte monocyte_IFN24.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt

java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in Fairfax_2014\
 --out monocyte_LPS2\
 --inexp monocyte_LPS2.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte monocyte_LPS2.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt


 java -Xmx10g -Xms10g -XX:StringTableSize=319973\
 -jar ~/software/eqtl-mapping-pipeline-1.4.7-SNAPSHOT/eqtl-mapping-pipeline.jar\
 --mode mixupmapper\
 --in Fairfax_2014\
 --out B-cell_CD19\
 --inexp B-cell_CD19.tsv\
 --inexpplatform HT12v4\
 --inexpannot phenotype_metadata.tsv\
 --gte B-cell_CD19.geno_pheno_coupling.tsv\
 --testall\
 --snps snp_list.txt