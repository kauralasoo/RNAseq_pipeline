rule run_CrossMap:
    input:
        vcf = "GRCh37/_EGAZ00001235598_release_vcf_06092016_All_chr.BPWP10_13_12_15.vcf.gz",
        chain = "GRCh37_to_GRCh38.chain",
        ref_genome = "/gpfs/rocket/home/a72094/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
    output:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.vcf.gz",
    params:
        temp_vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.vcf"
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

rule postprocess_CrossMap:
    input:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.vcf.gz"
    output:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.post.vcf.gz"
    threads: 1
    resources:
        mem = 1000
    shell:
        """
        source activate py2.7
        module load bcftools-1.6
        python postprocessCrossmap.py --vcf {input.vcf} | bgzip > {output.vcf}
        """

rule sort_vcf:
    input:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.post.vcf.gz"
    output:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.vcf.gz"
    threads: 1
    resources:
        mem = 35000
    shell:
        """
        module load bcftools-1.6
        bcftools sort -m 20000M -o {output.vcf} -O z {input.vcf}
        """

rule filter_ref_allele:
    input:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.vcf.gz",
        fasta = "/gpfs/rocket/home/a72094/annotations/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
    output:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.ref.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        module load bcftools-1.6
        bcftools norm -c x -O z -f {input.fasta} {input.vcf} > {output.vcf}
        """

rule remove_multialleic:
    input:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.ref.vcf.gz"
    output:
        vcf = "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.ref.filtered.vcf.gz"
    threads: 1
    resources:
        mem = 2000
    shell:
        """
        module load bcftools-1.6
        bcftools norm -m+any {input.vcf} | bcftools view -m2 -M2 - | bcftools annotate --set-id +'%CHROM\_%POS' | bcftools norm -d both -O z > {output.vcf}
        """

rule keep_common:
    input:
        "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.ref.filtered.vcf.gz"
    output:
        "GRCh38/BLUEPRINT_06092016_GRCh38.sorted.MAF05.vcf.gz"
    threads: 1
    resources:
        mem = 3000
    shell:
        """
        module load bcftools-1.6
        bcftools filter -i 'MAF[0] >= 0.05' -O z {input} > {output}
        """
