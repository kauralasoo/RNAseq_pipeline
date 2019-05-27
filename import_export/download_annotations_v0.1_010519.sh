#GENCODE version: 30
#Ensembl version: 96

#Download GENCODE annotations
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_30/gencode.v30.annotation.gtf.gz
mv gencode.v30.annotation.gtf.gz raw

#Rename chromosome names to Ensembl
zcat raw/gencode.v30.annotation.gtf.gz | sed 's/chrM/chrMT/;s/chr//' | gzip > gencode.v30.annotation.no_chr.gtf

#Download Ensembl gene annotations
wget ftp://ftp.ensembl.org/pub/release-96/gtf/homo_sapiens/Homo_sapiens.GRCh38.96.gtf.gz
gunzip Homo_sapiens.GRCh38.96.gtf.gz

#Download the reference genome
wget ftp://ftp.ensembl.org/pub/release-96/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
gunzip Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz

#Index the referece genome
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa

#Extract exons and splice sites to build hisat2 index (HISAT2 version 2.1.0)
mkdir hisat2_index_v96
hisat2_extract_splice_sites.py Homo_sapiens.GRCh38.96.gtf > Homo_sapiens.GRCh38.96.splice_sites.txt
hisat2_extract_exons.py Homo_sapiens.GRCh38.96.gtf > Homo_sapiens.GRCh38.96.exons.txt
hisat2-build --ss Homo_sapiens.GRCh38.96.splice_sites.txt --exon Homo_sapiens.GRCh38.96.exons.txt Homo_sapiens.GRCh38.dna.primary_assembly.fa hisat2_index_v96/Homo_sapiens.GRCh38.dna.primary_assembly

#Extract intron information for Leafcutter

#Download txrevise annotations
wget https://zenodo.org/record/3232932/files/Homo_sapiens.GRCh38.96.version_1.tar.gz
tar xzfv Homo_sapiens.GRCh38.96.version_1.tar.gz

