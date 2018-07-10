CHROMS = ["1","10","11","12","13","14","15","16","17","18","19","2","20","21","22","3","4","5","6","7","8","9"]

rule rename_variants:
    input:
        gen = "/gpfs/hpchome/a72094/rocket/datasets/TwinsUK/genotypes/gen/chr{chrom}_Eurobats_Public.gen"
    output:
        gen = "processed/gen/chr_{chrom}.gen"
    shell:
        "../../genotype_scripts/gen_replace_variant_id.py --gen {input.gen} --chr {wildcards.chrom} > {output.gen}"

rule make_vcf:
    input:
        gen = "processed/gen/chr_{chrom}.gen"
        sample = "/gpfs/hpchome/a72094/rocket/datasets/TwinsUK/genotypes/sample_lists/chr{chrom}_Eurobats_Public.sample"
    output:
        vcf = "processed/vcf/chr_{chrom}.vcf"
    shell:
        """
        module load bcftools-1.8
        bcftools convert --gensample2vc {input.gen},{input.sample} -Oz -o {output.vcf}
        """

rule merge_vcfs:
    input:
        expand("processed/vcf/chr_{chrom}.vcf", chrom = CHROM)
    output:
        "processed/TwinsUK_GRCh37_1KG_Phase1.vcf.gz"
    shell:
        """
        module load bcftools-1.8
        bcftools concat {input} -Oz -o {output}
        """
