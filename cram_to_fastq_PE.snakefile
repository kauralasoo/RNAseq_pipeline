import uuid
import os

#Convert BAM files to fastq
rule cram_to_fastq:
    input:
        lambda wildcards: config["samples"][wildcards.sample]
    output:
        fq1 = "processed/{dataset}/fastq/{sample}/{sample}_1.fastq.gz",
        fq2 = "processed/{dataset}/fastq/{sample}/{sample}_2.fastq.gz",
    params:
        reference_annotation = "/gpfs/hpchome/a72094/rocket/annotations/GRCh37/hs37d5.fa",
        local_tmp = "/tmp/kerimov_" + uuid.uuid4().hex + "/"
    resources:
        mem = 3000
    threads: 1
    shell:
        """
        module load samtools-1.6
        mkdir {params.local_tmp}
        rsync -aP --bwlimit=10000 {input[0]} {params.local_tmp}/{wildcards.sample}.cram
        rsync -aP --bwlimit=10000 {params.reference_annotation} {params.local_tmp}/reference_annotation_tmp
        samtools collate --reference {params.local_tmp}/reference_annotation_tmp {params.local_tmp}/{wildcards.sample}.cram {params.local_tmp}/{wildcards.sample}.collated
        samtools fastq -F 2816 -c 6 -1 {params.local_tmp}/{wildcards.sample}_1.fastq.gz -2 {params.local_tmp}/{wildcards.sample}_2.fastq.gz {params.local_tmp}/{wildcards.sample}.collated.bam
        rsync -aP --bwlimit=10000 {params.local_tmp}/{wildcards.sample}_1.fastq.gz {output.fq1}
        rsync -aP --bwlimit=10000 {params.local_tmp}/{wildcards.sample}_2.fastq.gz {output.fq2}
        rm -r {params.local_tmp}
        """

#Make sure that all final output files get created
rule make_all:
    input:
        expand("processed/{{dataset}}/fastq/{sample}/{sample}_1.fastq.gz", sample=config["samples"])
    output:
        "processed/{dataset}/out.txt"
    resources:
        mem = 100
    threads: 1
    shell:
        "echo 'Done' > {output}"

# samtools collate --reference /tmp/kerimov/hs37d5.fa processed/SensoryNeurons/PART_3/EGAR00001434763_18676_4#21.cram processed/SensoryNeurons/PART_3/EGAR00001434763_18676_4#21.collated