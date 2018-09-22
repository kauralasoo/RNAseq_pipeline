rule crossmap_genotypes:
    input:
        vcf = expand("processed/{{study}}/genotypes/GRCh38/chr{chromosome}.dose.vcf.gz", chromosome = config["chromosomes"])
    output:
        "processed/{study}/out.txt"
    threads: 1
    resources:
        mem = 1000
    shell:
        "echo 'Done! > {output}'" 

rule run_CrossMap:
    input:
        vcf = "processed/{study}/genotypes/GRCh37/chr{chromosome}.dose.vcf.gz",
    output:
        vcf = "processed/{study}/genotypes/GRCh38/chr{chromosome}.dose.vcf.gz",
    params:
        temp_vcf = "processed/{study}/genotypes/GRCh38/chr{chromosome}.dose.vcf"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        source activate py2.7
        module load samtools-1.6
        CrossMap.py vcf {config[chain_file]} {input.vcf} {config[ref_genome]} {params.temp_vcf}
        bgzip {params.temp_vcf}
        """