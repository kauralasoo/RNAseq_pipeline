CHROMS = ["1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9","X","Y"]

rule make_all:
    input:
        expand("/gpfs/hpchome/a72094/hpc/datasets/controlled_access/Garieri_2017/1000G/Garieri_1000G.chr{chromosome}.vcf.gz", chromosome = CHROMS)
    output:
        "out.txt"
    threads: 1
    resources:
        mem = 1000
    shell:
        "echo 'Done!' > {output}" 
        
rule extract_genotypes:
    input:
        vcf = "/gpfs/hpchome/a72094/rocket/datasets/1000G/GRCh38/ALL.chr{chromosome}_GRCh38.genotypes.20170504.vcf.gz",
        samples = "/gpfs/hpchome/a72094/hpc/projects/RNAseq_pipeline/metadata/Garieri_2017/Garieri_genotype_names.txt"
    output:
        vcf = "/gpfs/hpchome/a72094/hpc/datasets/controlled_access/Garieri_2017/1000G/Garieri_1000G.chr{chromosome}.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        module load bcftools-1.9
        bcftools view -S {input.samples} --force-samples {input.vcf} -Oz -o {output.vcf}
        """

