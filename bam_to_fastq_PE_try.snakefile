import uuid
import os

#Convert BAM files to fastq
rule bam_to_fastq:
	input:
		"Fat/{sample}.bam"
	output:
		fq1 = "processed/{dataset}/fastq/{sample}_1.fastq.gz",
		fq2 = "processed/{dataset}/fastq/{sample}_2.fastq.gz",
	resources:
		mem = 3000
	threads: 1
	shell:
		"""
		module load samtools-1.6
		samtools collate {input} collatedFiles/{wildcards.sample}.collated
		samtools fastq -F 2816 -c 6 -1 {output.fq1} -2 {output.fq2} collatedFiles/{input}.collated.bam
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{dataset}}/fastq/{sample}_1.fastq.gz", sample=config["samples"])
	output:
		"processed/{dataset}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"
