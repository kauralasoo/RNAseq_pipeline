rule run_CrossMap:
    input:
        vcf = "/gpfs/hpchome/a72094/rocket/datasets/TwinsUK/genotypes/GRCh37/TwinsUK_GRCh37_1KG_Phase1.vcf.gz",
        chain = "../../GRCh37_to_GRCh38.chain",
        ref_genome = "/gpfs/rocket/home/a72094/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        vcf = "processed/TwinsUK_GRCh38.vcf.gz",
    params:
        temp_vcf = "processed/TwinsUK_GRCh38.vcf"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        source activate py2.7
        module load samtools-1.6
        CrossMap.py vcf {input.chain} {input.vcf} {input.ref_genome} {params.temp_vcf}
        bgzip {params.temp_vcf}
        """