rule map_qtls:
	input:
		expand("processed/{{study}}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz", annot_type = config["quant_methods"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/tab/{condition}.nominal.txt.gz", annot_type = config["quant_methods"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz", annot_type = config["quant_methods"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz.tbi", annot_type = config["quant_methods"], condition = config["conditions"]),
		expand("processed/{{study}}/qtltools/output/{annot_type}/final/{condition}.variant_information.txt.gz", annot_type = config["quant_methods"], condition = config["conditions"]),
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Compres and index input bed file
rule compress_bed:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed"
	output:
		bed = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz"),
		bed_index = protected("processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi")
	threads: 1
	resources:
		mem = 100
	shell:
		"""
		module load samtools-1.6
		bgzip {input.bed} && tabix -p bed {output.bed}
		"""

#Extract samples from vcf
rule extract_samples:
	input:
		samples = "processed/{study}/qtltools/input/{annot_type}/{condition}.sample_names.txt"
	output:
		vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
		vcf_index = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz.csi"
	threads: 1
	resources:
		mem = 100
	shell:
		"""
		module load bcftools-1.8
		bcftools view -S {input.samples} {config[vcf_file]} -Oz -o {output.vcf}
		bcftools index {output.vcf}
		"""

#Extract variant information from VCF
rule extract_variant_information:
	input:
		vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
	output:
		var_info = "processed/{study}/qtltools/output/{annot_type}/final/{condition}.variant_information.txt.gz"
	threads: 1
	resources:
		mem = 100
	run:
		shell("module load bcftools-1.8")
		if(config["is_imputed"] == False):
			shell("bcftools +fill-tags {input.vcf} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\tNA\n' | gzip > {output.var_info}")
		else:
			shell("bcftools +fill-tags {input.vcf} | bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%TYPE\t%AC\t%AN\t%MAF\tR2\n' | gzip > {output.var_info}")

#Perform PCA on the genotype and phenotype data
rule perform_pca:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi",
		vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz",
		vcf_index = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz.csi"
	output:
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates.txt"
	params:
		pheno_pca = "processed/{study}/qtltools/input/{annot_type}/{condition}.pheno",
		geno_pca = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.geno"
	threads: 1
	resources:
		mem = 2000
	shell:
		"""
		QTLtools pca --bed {input.bed} --center --scale --out {params.pheno_pca}
		QTLtools pca --vcf {input.vcf} --maf 0.05 --center --scale --distance 50000 --out {params.geno_pca}
		head -n 7 {params.pheno_pca}.pca > {output.covariates}
		set +o pipefail; tail -n+2 {params.geno_pca}.pca | head -n 6 >> {output.covariates}
		"""


#Run QTLtools in permutation mode
rule permutation_run:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi",
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates.txt",
		vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz"
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/batches/{condition}.permutation.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[cis_window]} --permute 10000 --grp-best"


#Merge all batches from QTLtools
rule merge_permutation_batches:
	input:
		expand("processed/{{study}}/qtltools/output/{{annot_type}}/batches/{{condition}}.permutation.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		"processed/{study}/qtltools/output/{annot_type}/{condition}.permuted.txt.gz"
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		module load samtools-1.6
		cat {input} | bgzip > {output}
		"""


#Run QTLtools in nominal mode
rule nominal_run:
	input:
		bed = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz",
		bed_index = "processed/{study}/qtltools/input/{annot_type}/{condition}.bed.gz.tbi",
		covariates = "processed/{study}/qtltools/input/{annot_type}/{condition}.covariates.txt",
		vcf = "processed/{study}/qtltools/input/{annot_type}/vcf/{condition}.vcf.gz"
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/nominal_batches/{condition}.nominal.batch.{batch}.{n_batches}.txt")
	params:
		chunk = "{batch} {n_batches}"
	threads: 1
	resources:
		mem = 5000
	shell:
		"QTLtools cis --vcf {input.vcf} --bed {input.bed} --cov {input.covariates} --chunk {params.chunk} --out {output} --window {config[nominal_cis_window]} --nominal 1"

#Merge all batches from QTLtools
rule merge_nominal_batches:
	input:
		expand("processed/{{study}}/qtltools/output/{{annot_type}}/nominal_batches/{{condition}}.nominal.batch.{batch}.{n_batches}.txt", 
			batch=[i for i in range(1, config["n_batches"] + 1)],
			n_batches = config["n_batches"])
	output:
		temp("processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz")
	resources:
		mem = 100
	threads: 1
	shell:
		"""
		module load samtools-1.6
		cat {input} | bgzip > {output}
		"""

#Replace tabs
rule replace_space_tabs:
	input:
		"processed/{study}/qtltools/output/{annot_type}/{condition}.nominal.txt.gz"
	output:
		"processed/{study}/qtltools/output/{annot_type}/tab/{condition}.nominal.txt.gz"
	resources:
		mem = 1000
	threads: 2
	shell:
		"""
		gzip -dc {input} | awk -v OFS='\\t' '{{$1=$1; print $0}}' | gzip > {output}
		"""

#Add SNP coordinates to QTLTools output file
rule sort_qtltools_output:
	input:
		"processed/{study}/qtltools/output/{annot_type}/tab/{condition}.nominal.txt.gz"
	output:
		protected("processed/{study}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz")
	resources:
		mem = 12000
	threads: 10
	shell:
		"""
		module load samtools-1.6
		gzip -dc {input} | LANG=C sort -k9,9 -k10,10n -k11,11n -S11G --parallel=8 | bgzip > {output}
		"""

#Tabix-index QTLtools output files
rule index_qtltools_output:
	input:
		"processed/{study}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz"
	output:
		"processed/{study}/qtltools/output/{annot_type}/final/{condition}.nominal.sorted.txt.gz.tbi"
	resources:
		mem = 1000
	threads: 1
	shell:
		"""
		module load samtools-1.6
		tabix -s9 -b10 -e11 -f {input}
		"""
