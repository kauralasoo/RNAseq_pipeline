import uuid
import os

#Run MBV on all samples
rule run_samtools_index:
    input:
        bam = "processed/{study}/hisat2/{sample}.bam"
    output:
        "processed/{study}/hisat2/{sample}.bam.bai"
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
        samtools index {params.local_tmp}/{wildcards.sample}.bam
        rsync -aP --bwlimit=10000 {params.local_tmp}/{wildcards.sample}.bam.bai {output}
        rm -r {params.local_tmp}
        """

#Make sure that all final output files get created
rule make_all:
    input:
        expand("processed/{{study}}/hisat2/{sample}.bam.bai", sample=config["samples"])
    output:
        "processed/{study}/out.txt"
    resources:
        mem = 100
    threads: 1
    shell:
        "echo 'Done' > {output}"
