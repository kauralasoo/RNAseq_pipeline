import uuid
import os

#Align reads to the reference genome using HISAT2
rule hisat2_align:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		bam = "processed/{study}/hisat2/{sample}.bam",
		ss = "processed/{study}/hisat2_ss/{sample}.splice_sites.txt"
	params:
		local_tmp = "/tmp/kerimov_" + uuid.uuid4().hex + "/"
	resources:
		mem = 8000
	threads: 8
	shell:
		"""
		module load samtools-1.6
		mkdir {params.local_tmp}
		rsync -aP --bwlimit=10000 {input[0]} {params.local_tmp}/{wildcards.sample}_1.fq.gz
		rsync -aP --bwlimit=10000 {input[1]} {params.local_tmp}/{wildcards.sample}_2.fq.gz
		hisat2 -p {threads} -x {config[hisat2_index]} {config[hisat2_flags]} --novel-splicesite-outfile {output.ss} -1 {params.local_tmp}/{wildcards.sample}_1.fq.gz -2 {params.local_tmp}/{wildcards.sample}_2.fq.gz | samtools view -Sb > {params.local_tmp}/{wildcards.sample}.bam
		samtools sort -m 1000M -o {params.local_tmp}/{wildcards.sample}.sorted.bam -O BAM --threads 6 {params.local_tmp}/{wildcards.sample}.bam
		rsync -aP --bwlimit=10000 {params.local_tmp}/{wildcards.sample}.sorted.bam {output.bam}
		rm -r {params.local_tmp}
		"""

#Sort bam files by name and quantify gene expression using featureCounts
rule quantify_featureCounts:
	input:
		bam = "processed/{study}/hisat2/{sample}.bam"
	output:
		counts = "processed/{study}/featureCounts/{sample}.featureCounts.txt",
		summary = "processed/{study}/featureCounts/{sample}.featureCounts.txt.summary"
	params:
		local_tmp = "/tmp/kerimov_" + uuid.uuid4().hex + "/"
	threads: 6
	resources:
		mem = 10000
	run:
		shell("mkdir {params.local_tmp}")
		shell("rsync -aP --bwlimit=10000 {input.bam} {params.local_tmp}/{wildcards.sample}.bam")
		shell("module load samtools-1.6 && samtools sort -n -m 1000M -o {params.local_tmp}/{wildcards.sample}.sorted.bam -O BAM --threads 5 {params.local_tmp}/{wildcards.sample}.bam")
		if(config["strandedness"] == "Stranded"):
			shell("featureCounts -s2 -p -C -D 5000 -d 50 --donotsort -a {config[ensembl_gtf]} -o {output.counts} {params.local_tmp}/{wildcards.sample}.sorted.bam")
		else:
			shell("featureCounts -s0 -p -C -D 5000 -d 50 --donotsort -a {config[ensembl_gtf]} -o {output.counts} {params.local_tmp}/{wildcards.sample}.sorted.bam")
		shell("rm -r {params.local_tmp}")
		
#Merge featureCounts results
rule merge_featureCounts:
	input:
		expand("processed/{{study}}/featureCounts/{sample}.featureCounts.txt", sample=config["samples"])
	output:
		"processed/{study}/matrices/gene_expression_featureCounts.txt"
	params:
		sample_ids = ','.join(config["samples"]),
		dir = "processed/{study}/featureCounts/"
	threads: 1
	resources:
		mem = 12000
	shell:
		"""
		module load R-3.4.1
		Rscript scripts/merge_featureCounts.R -s {params.sample_ids} -d {params.dir} -o {output}
		"""

#Run two leafcutter steps in one command
rule leafcutter_bam_to_junc:
	input:
		"processed/{study}/hisat2/{sample}.bam"
	output:
		"processed/{study}/leafcutter/junc/{sample}.junc"
	params:
		local_tmp = "/tmp/kerimov_" + uuid.uuid4().hex + "/"
	threads: 1
	resources:
		mem = 1000
	shell:
		"""
		module load samtools-1.6
		module load perl-5.22.0
		mkdir {params.local_tmp}
		rsync -aP --bwlimit=10000 {input} {params.local_tmp}/{wildcards.sample}.bam
		samtools view {params.local_tmp}/{wildcards.sample}.bam | python {config[leafcutter_root]}/scripts/filter_cs.py | {config[leafcutter_root]}/scripts/sam2bed.pl --use-RNA-strand - {params.local_tmp}/{wildcards.sample}.bed
		{config[leafcutter_root]}/scripts/bed2junc.pl {params.local_tmp}/{wildcards.sample}.bed {params.local_tmp}/{wildcards.sample}.junc
		cp {params.local_tmp}/{wildcards.sample}.junc {output}
		rm -r {params.local_tmp}
		"""

#Convert gff3 files generated by reviseAnnotations into fasta sequence
rule convert_gff3_to_fasta:
	input:
		"processed/annotations/gff/{annotation}.gff3"
	output:
		"processed/annotations/fasta/{annotation}.fa"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		module load cufflinks2.2
		gffread -w {output} -g {config[reference_genome]} {input}
		"""

#Build salmon indexes for fasta files
rule construct_salmon_index:
	input:
		"processed/annotations/fasta/{annotation}.fa"
	output:
		"processed/annotations/salmon_index/{annotation}"
	resources:
		mem = 10000
	threads: 1
	shell:
		"""
		salmon -no-version-check index -t {input} -i {output}
		"""

#Quantify gene expression using full Ensembl annotations
rule quant_salmon:
	input:
		lambda wildcards: config["samples"][wildcards.sample],
		"processed/annotations/salmon_index/{annotation}"
	output:
		protected("processed/{study}/salmon/{annotation}/{sample}/quant.sf")
	params:
		out_prefix = "processed/{study}/salmon/{annotation}/{sample}",
		local_tmp = "/tmp/" + uuid.uuid4().hex + "/"
	resources:
		mem = 10000
	threads: 8	
	shell:
		"""
		mkdir {params.local_tmp}
		mkdir {params.local_tmp}/salmon_index
		rsync -aP --bwlimit=10000 {input[2]}/* {params.local_tmp}/salmon_index/
		rsync -aP --bwlimit=10000 {input[0]} {params.local_tmp}/{wildcards.sample}_1.fq.gz
		rsync -aP --bwlimit=10000 {input[1]} {params.local_tmp}/{wildcards.sample}_2.fq.gz
		salmon --no-version-check quant --seqBias --useVBOpt --gcBias --libType {config[libType]} --index {params.local_tmp}/salmon_index -1 {params.local_tmp}/{wildcards.sample}_1.fq.gz -2 {params.local_tmp}/{wildcards.sample}_2.fq.gz -p {threads} -o {params.out_prefix}
		rm -r {params.local_tmp}
		"""

#Merge Salmon results
rule merge_salmon:
	input:
		expand("processed/{{study}}/salmon/{{annotation}}/{sample}/quant.sf", sample=config["samples"])
	output:
		"processed/{study}/matrices/{annotation}.salmon_txrevise.rds"
	params:
		sample_ids = ','.join(config["samples"]),
		dir = "processed/{study}/salmon/{annotation}"
	threads: 1
	resources:
		mem = 12000
	shell:
		"""
		module load R-3.4.1
		Rscript scripts/merge_Salmon.R -s {params.sample_ids} -d {params.dir} -o {output}
		"""

#Run MBV on all samples
rule run_qtltools_mbv:
	input:
		bam = "processed/{study}/hisat2/{sample}.bam"
	output:
		"processed/{study}/mbv/{sample}.mbv_output.txt"
	params:
		local_tmp = "/tmp/a72094_" + uuid.uuid4().hex + "/"
	resources:
		mem = 10000
	threads: 8
	shell:
		"""
		module load samtools-1.6
		mkdir {params.local_tmp}
		rsync -aP --bwlimit=10000 {input.bam} {params.local_tmp}/{wildcards.sample}.bam
		samtools index {params.local_tmp}/{wildcards.sample}.bam
		QTLtools mbv --vcf {config[vcf_file]} --bam {params.local_tmp}/{wildcards.sample}.bam --out {output}
		rm -r {params.local_tmp}
		"""
		
#Make sure that all final output files get created
rule make_all:
	input:
		expand("processed/{{study}}/mbv/{sample}.mbv_output.txt", sample=config["samples"]),
		#expand("processed/{{study}}/bigwig/{sample}.str1.bw", sample=config["samples"]),
		#expand("processed/{{study}}/hisat2/{sample}.bam", sample=config["samples"]),
		"processed/{{study}}/matrices/gene_expression_featureCounts.txt",
		# expand("processed/{{study}}/leafcutter/junc/{sample}.junc", sample=config["samples"]),
		#expand("processed/{study}/salmon/{annotation}/{sample}/quant.sf", study = config["study"], annotation=config["annotations"], sample=config["samples"]),
		#expand("processed/{{study}}/featureCounts/{sample}.featureCounts.txt", sample=config["samples"]),
		#expand("processed/{study}/ASEcounts/{sample}.ASEcounts", study = config["study"], sample=config["samples"]),
		# expand("processed/{{study}}/matrices/{annotation}.salmon_txrevise.rds", annotation=config["annotations"]),
		#"processed/{study}/leafcutter/leafcutter_perind.counts.gz"
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"
