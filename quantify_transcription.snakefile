#Align reads to the reference genome using STAR
rule star_align:
	input:
		lambda wildcards: config["samples"][wildcards.sample]
	output:
		bam = "processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	params:
		prefix = "processed/{study}/STAR/{sample}/{sample}.",
		rg = 'ID:1 \"LB:1\tPL:Illumina\tSM:{sample}\tPU:1\"',
		tmp_fq1 = tempfile.mkstemp(),
		tmp_fq2 = tempfile.mkstemp(),
		star_temp = tempfile.mkdtemp()
	resources:
		mem = 42000
	threads: 8
	shell:
		"module load star-2.5.2 && "
		"cp {input[0]} {params.tmp_fq1} && "
		"cp {input[1]} {params.tmp_fq2} && "
		"STAR --runThreadN {threads} --outSAMtype BAM SortedByCoordinate --outWigType bedGraph "
		"--outWigNorm None --outWigStrand {config[outWigStrand]} --outSAMattrRGline {params.rg} "
		"--readFilesCommand zcat --genomeDir {config[star_index]} --limitBAMsortRAM 32000000000 "
		"--outSAMunmapped Within --outFileNamePrefix {params.prefix} --outTmpDir {params.star_tmp} "
		"--readFilesIn {params.tmp_fq1} {params.tmp_fq2}"

#Index sorted bams
rule index_bams:
	input:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
	resources:
		mem = 50
	threads: 1
	shell:
		"module load samtools-1.6 && "
		"samtools index {input}"

#Check genotype concordance between RNA-seq and VCF
rule check_genotype_concordance:
	input:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
	output:
		"processed/{study}/verifyBamID/{sample}.verifyBamID.bestSM"
	params:
		out_prefix = "processed/{study}/verifyBamID/{sample}.verifyBamID"
	resources:
		mem = 1500
	threads: 1
	shell:
		"verifyBamID --vcf {config[vcf_file]} --bam {input} --out {params.out_prefix} --best --ignoreRG"

#Sort bedgraph files before bigwig conversion
rule sort_bedgraph:
	input:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		bg = temp("processed/{study}/STAR/{sample}/{sample}.Signal.Unique.{strand}.sorted.bg")
	params:
		bg = "processed/{study}/STAR/{sample}/{sample}.Signal.Unique.{strand}.out.bg",
		bg_mult = "processed/{study}/STAR/{sample}/{sample}.Signal.UniqueMultiple.{strand}.out.bg"
	resources:
		mem = 1000
	threads: 1
	shell:
		"LC_COLLATE=C sort -k1,1 -k2,2n {params.bg} > {output.bg} && "
		"rm {params.bg} && rm {params.bg_mult}"


#Convert bedgraph to bigwig
rule bedgraph_to_bigwig:
	input:
		bg = "processed/{study}/STAR/{sample}/{sample}.Signal.Unique.{strand}.sorted.bg"
	output:
		bw = "processed/{study}/bigwig/{sample}.{strand}.bw",
	resources:
		mem = 1000
	threads: 1
	shell:
		"bedGraphToBigWig {input.bg} {config[chromosome_lengths]} {output.bw}"


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
		"module load cufflinks2.2 && "
		"gffread -w {output} -g {config[reference_genome]} {input}"

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
		"salmon -no-version-check index -t {input} -i {output}"

#Quantify gene expression using full Ensembl annotations
rule quant_salmon:
	input:
		lambda wildcards: config["samples"][wildcards.sample],
		"processed/annotations/salmon_index/{annotation}"
	output:
		protected("processed/{study}/salmon/{annotation}/{sample}/quant.sf")
	params:
		out_prefix = "processed/{study}/salmon/{annotation}/{sample}"
	resources:
		mem = 10000
	threads: 8	
	shell:
		"salmon --no-version-check quant --seqBias --gcBias --libType {config[libType]} "
		"--index {input[2]} -1 {input[0]} -2 {input[1]} -p {threads} "
		"-o {params.out_prefix}"

#Convert BAMs to bed for leafcutter
rule leafcutter_bam_to_bed:
	input:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		temp("processed/{study}/leafcutter/bed/{sample}.bed")
	threads: 3
	resources:
		mem = 1000
	shell:
		"module load samtools-1.6 && "
		"samtools view {input} | python {config[leafcutter_root]}/scripts/filter_cs.py | {config[leafcutter_root]}/scripts/sam2bed.pl --use-RNA-strand - {output}"

#Convert bed file to junctions
rule leadcutter_bed_to_junc:
	input:
		"processed/{study}/leafcutter/bed/{sample}.bed"
	output:
		temp("processed/{study}/leafcutter/junc/{sample}.junc")
	threads: 1
	resources:
		mem = 1000
	shell:
		"{config[leafcutter_root]}/scripts/bed2junc.pl {input} {output}"

#Cluster junctions with LeafCutter
rule leafcutter_cluster_junctions:
	input:
		expand("processed/{{study}}/leafcutter/junc/{sample}.junc", sample = config["samples"])
	output:
		"processed/{study}/leafcutter/leafcutter_perind.counts.gz"
	params:
		junc_files = "processed/{study}/leafcutter/junction_files.txt",
		out_prefix = "processed/{study}/leafcutter/"
	threads: 1
	resources:
		mem = 1000
	shell:
		"ls --color=never {params.out_prefix}/junc/*.junc | cat > {params.junc_files} && "
		"python {config[leafcutter_root]}/clustering/leafcutter_cluster.py -j {params.junc_files} -r {params.out_prefix} -m 50 -l 500000"

#Sort BAMs by name
rule sort_bam_by_name:
	input:
		"processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam"
	output:
		temp("processed/{study}/sorted_bam/{sample}.Aligned.sortedByName.out.bam")
	threads: 6
	resources:
		mem = 8000
	shell:
		"module load samtools-1.6 && "
		"samtools sort -n -m 1000M -o {output} -O BAM --threads 5 {input}"


#Quantify expression using featureCounts
rule quantify_featureCounts:
	input:
		bam = "processed/{study}/sorted_bam/{sample}.Aligned.sortedByName.out.bam"
	output:
		counts = "processed/{study}/featureCounts/{sample}.featureCounts.txt",
		summary = temp("processed/{study}/featureCounts/{sample}.featureCounts.txt.summary")
	threads: 1
	resources:
		mem = 1000
	shell:
		"featureCounts -p -C -D 5000 -d 50 --donotsort -a {config[ensembl_gtf]} -o {output.counts} {input.bam}"


#Quantify allele-specific expression
rule count_ASE:
	input:
		bam = "processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam",
		index = "processed/{study}/STAR/{sample}/{sample}.Aligned.sortedByCoord.out.bam.bai"
	output:
		counts = "processed/{study}/ASEcounts/{sample}.ASEcounts"
	resources:
		mem = 8000
	threads: 1
	params:
		gatk_command = "/software/java/bin/java -jar -Xmx6g ~/software/GenomeAnalysisTK.jar"
	shell:
		"{params.gatk_command} -T ASEReadCounter -R {config[fasta]} -I {input.bam} -o {output.counts} -sites {config[ase_vcf]} -U ALLOW_N_CIGAR_READS -dt NONE --minMappingQuality 10 -rf MateSameStrand"


#Make sure that all final output files get created
rule make_all:
	input:
		#expand("processed/{study}/verifyBamID/{sample}.verifyBamID.bestSM", study = config["study"], sample=config["samples"]),
		expand("processed/{study}/bigwig/{sample}.{strand}.bw", study = config["study"], sample=config["samples"], strand = config["bigwig_strandsq"]),
		expand("processed/{study}/salmon/{annotation}/{sample}/quant.sf", study = config["study"], annotation=config["annotations"], sample=config["samples"]),
		#expand("processed/{study}/featureCounts/{sample}.featureCounts.txt", study = config["study"], sample=config["samples"]),
		#expand("processed/{study}/ASEcounts/{sample}.ASEcounts", study = config["study"], sample=config["samples"]),
		#"processed/{study}/leafcutter/leafcutter_perind.counts.gz"
	output:
		"processed/{study}/out.txt"
	resources:
		mem = 100
	threads: 1
	shell:
		"echo 'Done' > {output}"





