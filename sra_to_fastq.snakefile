
#Convert BAM files to fastq
rule sra_to_fastq:
	input:
		"processed/{dataset}/SRR/{sample}.sra"
	output:
		fq1 = "processed/{dataset}/fastq/{sample}_1.fastq.gz"
	resources:
		mem = 1000
	params:
		outdir = "processed/{dataset}/fastq/"
	threads: 1
	shell:
		"""
		fastq-dump --split-files --gzip --skip-technical --readids --dumpbase -v -O {params.outdir} --clip {input}
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
