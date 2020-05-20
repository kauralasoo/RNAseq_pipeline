
# Add standard variant IDs
bcftools annotate BLUEPRINT_06092016_GRCh38.sorted.ref.filtered.vcf.gz --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o BLUEPRINT_06092016_GRCh38.annotated.vcf.gz

#Keep only unrelated individuals
bcftools filter BLUEPRINT_06092016_GRCh38.annotated.vcf.gz -i 'MAF[0] > 0.01' -Oz -o BLUEPRINT_06092016_GRCh38_filtered.vcf.gz
bcftools index BLUEPRINT_06092016_GRCh38_filtered.vcf.gz

#Extract variant information
module load bcftools-1.8
bcftools +fill-tags BLUEPRINT_06092016_GRCh38_filtered.vcf.gz | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\tNA\n' | gzip > BLUEPRINT_06092016_GRCh38.variant_information.txt.gz


### BLUEPRINT V2 genotypes ###
# Using CrossMap v0.41
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/10.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr10.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/11.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr11.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/12.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr12.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/13.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr13.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/14.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr14.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/15.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr15.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/16.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr16.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/17.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr17.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/18.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr18.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/19.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr19.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/1.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr1.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/20.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr20.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/21.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr21.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/22.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr22.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/2.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr2.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/3.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr3.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/4.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr4.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/5.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr5.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/6.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr6.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/7.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr7.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/8.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr8.vcf
CrossMap.py vcf GRCh37_to_GRCh38.chain GRCh37/9.BPWP10_23_05_17.vcf.gz ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa GRCh38/chr9.vcf

#Compress all vcfs
ls *.vcf | xargs -n1 bgzip -f

#Merge all VCFs together
bcftools concat chr1.vcf.gz chr2.vcf.gz chr3.vcf.gz chr4.vcf.gz chr5.vcf.gz chr6.vcf.gz chr7.vcf.gz chr8.vcf.gz chr9.vcf.gz chr10.vcf.gz chr11.vcf.gz chr12.vcf.gz chr13.vcf.gz chr14.vcf.gz chr15.vcf.gz chr16.vcf.gz chr17.vcf.gz chr18.vcf.gz chr19.vcf.gz chr20.vcf.gz chr21.vcf.gz chr22.vcf.gz -Oz -o BLUEPRINT_merged.vcf.gz

#Filter
bcftools annotate .vcf.gz --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o BLUEPRINT_06092016_GRCh38.annotated.vcf.gz

#Filter
bcftools +fill-tags BLUEPRINT_merged.vcf.gz | bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' | bcftools norm -d all -Oz -o BLUEPRINT_filtered.vcf.gz

#Keep only relevant fields for eQTL analysis
bcftools annotate -x ^INFO/AR2,^FORMAT/GT,^FORMAT/DS BLUEPRINT_filtered.vcf.gz -Oz -o BLUEPRINT_pruned.vcf.gz

#Run checkref
bcftools norm --check-ref x -f ~/hpcproject/genomic_references/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa BLUEPRINT_pruned.vcf.gz -Oz -o BLUEPRINT_checkref.vcf.gz

#Set correct variant id
bcftools annotate BLUEPRINT_checkref.vcf.gz --set-id 'chr%CHROM\_%POS\_%REF\_%FIRST_ALT' -Oz -o BLUEPRINT.MAF001.vcf.gz

#Sort the vcf file by position
bcftools sort BLUEPRINT.MAF001.vcf.gz -Oz -o BLUEPRINT.MAF001.sorted.vcf.gz
