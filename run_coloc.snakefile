rule run_all:
	input:
		expand("results/coloc/{qtl_set}.{gwas_trait}.txt", qtl_set = config["qtl_sets"], gwas_trait = config["gwas_traits"])
	output:
		"results/coloc/coloc_out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done!' > {output}"

#Run coloc accross all traits
rule run_coloc:
	input:
		lead_vars = "results/qtl_summary_stats/{qtl_set}.permuted.txt.gz",
		summary_stats = "results/qtl_summary_stats/{qtl_set}.permuted.txt.gz",
		var_info = "results/qtl_summary_stats/{qtl_set}.variant_information.txt.gz",
	output:
		"results/coloc/{qtl_set}.{gwas_trait}.txt"
	resources:
		mem = 12000
	threads: 1
	shell:
		"""
		module load R/3.5.1
		Rscript ../colocWrapper/examples/qtlmap_run_coloc.R	-l {input.lead_vars} -w {config[coloc_window]} --gwas {wildcards.gwas_trait} -d {config[gwas_dir]} -o {output} --qtl {input.summary_stats} --qtlvarinfo {input.var_info} --gwaslist {config[gwas_list]}
		"""
