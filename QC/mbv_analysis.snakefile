import uuid
import os

#Run MBV on all samples
rule run_qtltools_mbv:
    input:
        bam = "processed/{dataset}/hisat2/{sample}.bam"
    output:
        "processed/{dataset}/mbv/{sample}.mbv_output.txt"
    resources:
        mem = 3000
    params:
        local_tmp = "/tmp/kerimov_" + uuid.uuid4().hex + "/"
    threads: 1
    shell:
        """
        module load samtools-1.6
        mkdir {params.local_tmp}
        rsync -aP --bwlimit=10000 {input.bam} {params.local_tmp}/{wildcards.sample}.bam
        rsync -aP --bwlimit=10000 {input.bam}.bai {params.local_tmp}/{wildcards.sample}.bam.bai
        QTLtools mbv --vcf {config[vcf_file]} --bam {params.local_tmp}/{wildcards.sample}.bam --out {params.local_tmp}/{wildcards.sample}.mbv_output.txt
        rsync -aP --bwlimit=10000 {params.local_tmp}/{wildcards.sample}.mbv_output.txt {output}
        rm -r {params.local_tmp}
        """

#Make sure that all final output files get created
rule make_all:
    input:
        expand("processed/{{dataset}}/mbv/{sample}.mbv_output.txt", sample=config["samples"])
    output:
        "processed/{dataset}/out.txt"
    resources:
        mem = 100
    threads: 1
    shell:
        "echo 'Done' > {output}"
