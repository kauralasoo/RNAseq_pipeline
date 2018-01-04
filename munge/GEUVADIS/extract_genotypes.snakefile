CHROMS = ["1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y"]

rule make_all:
    input:
        expand("/gpfs/hpchome/a72094/rocket/projects/GEUVADIS/genotypes/GEUVADIS_445_samples.chr{chromosome}.vcf.gz", chromosome = CHROMS)
    output:
        "out.txt"
    threads: 1
    resources:
        mem = 1000
    shell:
        "echo 'Done! > {output}'" 
        
rule extract_genotypes:
    input:
        vcf = "/gpfs/hpchome/a72094/rocket/datasets/1000G/GRCh38/ALL.chr{chromosome}_GRCh38.genotypes.20170504.vcf.gz",
        samples = "../../metadata/GEUVADIS/GEUVADIS_genotype_list.txt"
    output:
        vcf = "/gpfs/hpchome/a72094/rocket/projects/GEUVADIS/genotypes/GEUVADIS_445_samples.chr{chromosome}.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        module load bcftools-1.6
        bcftools view -S {input.samples}  --force-samples {input.vcf} | bcftools filter -i 'MAF[0] >= 0.05' -O z - > {output.vcf}
        """

#bcftools concat GEUVADIS_445_samples.chr1.vcf.gz GEUVADIS_445_samples.chr2.vcf.gz GEUVADIS_445_samples.chr10.vcf.gz GEUVADIS_445_samples.chr11.vcf.gz GEUVADIS_445_samples.chr12.vcf.gz GEUVADIS_445_samples.chr13.vcf.gz GEUVADIS_445_samples.chr14.vcf.gz GEUVADIS_445_samples.chr15.vcf.gz GEUVADIS_445_samples.chr16.vcf.gz GEUVADIS_445_samples.chr17.vcf.gz GEUVADIS_445_samples.chr18.vcf.gz GEUVADIS_445_samples.chr19.vcf.gz GEUVADIS_445_samples.chr20.vcf.gz GEUVADIS_445_samples.chr21.vcf.gz GEUVADIS_445_samples.chr22.vcf.gz GEUVADIS_445_samples.chr3.vcf.gz GEUVADIS_445_samples.chr4.vcf.gz GEUVADIS_445_samples.chr5.vcf.gz GEUVADIS_445_samples.chr6.vcf.gz GEUVADIS_445_samples.chr7.vcf.gz GEUVADIS_445_samples.chr8.vcf.gz GEUVADIS_445_samples.chr9.vcf.gz -O z -o GEUVADIS_445_samples.merged.vcf.gz