
#Convert BAM files to fastq
rule bam_to_fastq:
	input:
		"processed/{dataset}/dummy_input.txt"
	output:
		fq1 = "processed/{dataset}/fastq/{sample}_1.fastq.gz",
		fq2 = "processed/{dataset}/fastq/{sample}_2.fastq.gz",
	resources:
		mem = 200
	params:
		input_string = lambda wildcards: config["samples"][wildcards.sample]
	threads: 1
	shell:
		"""
		echo -e {params.input_string} | python scripts/merge_fastq.py --indir processed/{dataset}/raw_fastq/ --outdir processed/{dataset}/fastq/ --suffix _1.fastq.gz
		echo -e {params.input_string} | python scripts/merge_fastq.py --indir processed/{dataset}/raw_fastq/ --outdir processed/{dataset}/fastq/ --suffix _2.fastq.gz
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
