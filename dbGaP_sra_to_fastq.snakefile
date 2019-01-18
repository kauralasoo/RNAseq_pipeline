
#Convert BAM files to fastq
rule sra_to_fastq:
	input:
		"sra/{sample}.sra"
	output:
		fq1 = "fastq/{sample}_1.fastq.gz"
	resources:
		mem = 1000
	params:
		outdir = "fastq/"
	threads: 1
	shell:
		"""
		fastq-dump --split-files --gzip --skip-technical --readids --dumpbase -v -O {params.outdir} --clip {input}
		"""

#Make sure that all final output files get created
rule make_all:
	input:
		expand("fastq/{sample}_1.fastq.gz", sample=config["srr_samples"])
	output:
		"out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"
