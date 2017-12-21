import uuid
import os

#Convert BAM files to fastq
rule bam_to_fastq:
	input:
		"processed/{dataset}/input_bam/{sample}.bam"
	output:
		fq1 = "processed/{dataset}/fastq/{sample}_1.fastq.gz",
		fq2 = "processed/{dataset}/fastq/{sample}_2.fastq.gz",
	params:
		local_tmp = "/tmp/" + uuid.uuid4().hex + "/"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		module load samtools-1.6
		mkdir {params.local_tmp}
		cp {input} {params.local_tmp}/{sample}.bam
		samtools {params.local_tmp}/{sample}.bam {params.local_tmp}/{sample}.collated
		samtools fastq -F 2816 -c 6 -1 {params.local_tmp}/{sample}_1.fastq.gz -2 {params.local_tmp}/{sample}_2.fastq.gz {params.local_tmp}/{sample}.collated.bam
		cp {params.local_tmp}/{sample}_1.fastq.gz {output.fq1}
		cp {params.local_tmp}/{sample}_2.fastq.gz {output.fq2}
		rm -r {params.local_tmp}
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
