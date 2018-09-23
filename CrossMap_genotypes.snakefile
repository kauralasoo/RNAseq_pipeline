rule crossmap_genotypes:
    input:
        vcf = expand("processed/{{study}}/genotypes/GRCh38/chr{chromosome}.dose.vcf.gz", chromosome = config["chromosomes"])
    output:
        "processed/{study}/out.txt"
    threads: 1
    resources:
        mem = 1000
    shell:
        "echo 'Done!' > {output}" 

rule run_CrossMap:
    input:
        vcf = "processed/{study}/genotypes/GRCh37/chr{chromosome}.dose.vcf.gz",
    output:
        vcf = temp("processed/{study}/genotypes/temp/chr{chromosome}.dose.vcf.gz"),
    params:
        temp_vcf = "processed/{study}/genotypes/temp/chr{chromosome}.dose.vcf"
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

rule filter_R2:
    input:
        vcf = "processed/{study}/genotypes/temp/chr{chromosome}.dose.vcf.gz"
    output:
        vcf = temp("processed/{study}/genotypes/filtered/chr{chromosome}.dose.vcf.gz")
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        module load bcftools-1.9
        bcftools filter -i 'INFO/R2 > {config[r2_thresh]}' {input.vcf} -Oz -o {output.vcf}
        """

rule sort_vcf:
    input:
        vcf = "processed/{study}/genotypes/filtered/chr{chromosome}.dose.vcf.gz"
    output:
        vcf = temp("processed/{study}/genotypes/sorted/chr{chromosome}.dose.vcf.gz"),
        csi = temp("processed/{study}/genotypes/sorted/chr{chromosome}.dose.vcf.gz.csi")
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        module load bcftools-1.9
        bcftools sort {input.vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """

rule keep_correct_chromosome:
    input:
        vcf = "processed/{study}/genotypes/sorted/chr{chromosome}.dose.vcf.gz",
        csi = "processed/{study}/genotypes/sorted/chr{chromosome}.dose.vcf.gz.csi"
    output:
        vcf = "processed/{study}/genotypes/GRCh38/chr{chromosome}.dose.vcf.gz",
        csi = "processed/{study}/genotypes/GRCh38/chr{chromosome}.dose.vcf.gz.csi"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        module load bcftools-1.9
        bcftools view -r {wildcards.chromosome} {input.vcf} -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """