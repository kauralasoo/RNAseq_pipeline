#Prepare genotypes
rule prepare_genotypes:
	input:
		bfile = expand("{file}.bed", file = config["raw_genotypes"])
	output:
		vcf = "processed/{study}.fixed.vcf.gz"
	threads: 4
	params:
		bfile = expand("{file}", file = config["raw_genotypes"])
	resources:
		mem = 4000
	shell:
		"""
		module load plink-1.9.0
		module load bcftools-1.8
		
		#Update to correct build
		update_build.sh {params.bfile} {config[strand_file]} processed/{wildcards.study}

		#Fix the ref and alt alleles
		plink -bfile processed/{wildcards.study} --reference-allele {config[refalt_file]} --make-bed --out processed/{wildcards.study}.RefAlt

		#Convert to VCF
		plink --bfile processed/{wildcards.study}.RefAlt --recode vcf-iid --output-chr M --out processed/{wildcards.study}.RefAlt

		#Fix ref and alt alleles using the dbSNP file
		bcftools view -t ^M processed/{wildcards.study}.RefAlt.vcf | bcftools +fixref - -Oz -o processed/{wildcards.study}.RefAlt.dbSNP.vcf.gz -- -f {config[ref_genome]} -i {config[dbSNP]}

		#Sort the vcf file
		bcftools sort processed/{wildcards.study}.RefAlt.dbSNP.vcf.gz -Oz -o processed/{wildcards.study}.RefAlt.dbSNP.sorted.vcf.gz

		#Remove variants that do not match reference, duplicates, and multi-alleics
		bcftools norm -cs -f {config[ref_genome]} processed/{wildcards.study}.RefAlt.dbSNP.sorted.vcf.gz | bcftools norm -d all - | bcftools norm -m+any - | bcftools view -m2 -M2 - -Oz -o processed/{wildcards.study}.RefAlt.dbSNP.sorted.filtered.vcf.gz

		#Calculate INFO tags
		bcftools +fill-tags processed/{wildcards.study}.RefAlt.dbSNP.sorted.filtered.vcf.gz | bcftools annotate --set-id +'%CHROM\_%POS' - -Oz -o processed/{wildcards.study}.RefAlt.dbSNP.sorted.filtered.tags.vcf.gz

		#Filter rare (MAF<0.01) and non-HWE varaints and those with abnormal reference alleles
		bcftools filter -i 'INFO/HWE > 1e-6 & F_MISSING < 0.05 & MAF[0] > 0.01' processed/{wildcards.study}.RefAlt.dbSNP.sorted.filtered.tags.vcf.gz -Ou | bcftools filter -e 'REF="N" | REF="I" | REF="D"' - -Oz -o {output}

		#Index the VCF file
		bcftools index {output}
		"""

rule extract_chromosome:
	input:
		vcf = "processed/{study}.fixed.vcf.gz"
	output:
		vcf = "processed/{study}/chr_{chr}.vcf.gz"
	threads: 1
	resources:
		mem = 1000
	shell:
		"""
		module load bcftools-1.8
		bcftools view -r {wildcards.chr} {input.vcf} -Oz o {output.vcf}
		"""

rule make_all:
	input:
		vcf = expand("processed/{{study}}/chr_{chr}.vcf.gz", chr = config["chromosomes"])
	output:
		"processed/out.txt"
	threads: 1
	resources:
		mem = 100
	shell:
		"""
		echo 'Done' > {output}
		"""
