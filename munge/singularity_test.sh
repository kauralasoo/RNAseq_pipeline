sudo singularity build --sandbox ubuntu shub://singularityhub/ubuntu
sudo singularity shell ubuntu/

apt-get update
apt-get install wget bzip2 -y

wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p /software/minconda3
rm Miniconda3-latest-Linux-x86_64.sh
PATH="/software/minconda3/bin:$PATH"

conda update -y conda
conda config --add channels defaults
conda config --add channels conda-forge
conda config --add channels bioconda

conda install -y salmon
conda install -y hisat2
conda install -y samtools

wget https://qtltools.github.io/qtltools/binaries/QTLtools_1.1_Ubuntu16.04_x86_64.tar.gz


export CPPFLAGS='-I/software/minconda3/include'
export LDFLAGS='-L/software/minconda3/lib'