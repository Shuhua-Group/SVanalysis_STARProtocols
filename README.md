# SVanalysis_STARProtocols
A protocol for applying low-coverage whole-genome sequencing data in structural variation studies


For complete details on the use and execution of this protocol, please refer to Liu et al. (2023)[1]


Reference


[1] Liu, Q., Yang, K., Xie, B., Gao, Y., Xu, S., and Lu, Y. (2023). Mapping Structural Variations in Haemaphysalis longicornis and Rhipicephalus microplus Reveals Vector–Pathogen Adaptation. iScience 26, 106398

## How to install docker image from tar

A docker image file can be downloaded from Google Drive (https://drive.google.com/drive/folders/10Ajzc8RfF9xMoRMLYuKnzSvGCjjO-wvP?usp=sharing)
```
#! 1- we have export a docker image file and upload to Google Drive
docker run -itd --name export svanalysis_starprotocols:1.0.0
docker export ddd44a6aa660 > svanalysis_starprotocols.tar
7z a docker_image_svanalysis_starprotocols.7z svanalysis_starprotocols.tar

#! 2- user can download and import our image
7z x docker_image_svanalysis_starprotocols.7z
docker import - svanalysis_starprotocols < svanalysis_starprotocols.tar
docker images

#! then you can see a image named "svanalysis_starprotocols"
```

## How to install docker image from docker file

1. go to the directory including the dockerfile
2. build image: docker build -t svanalysis_starprotocols:1.0.0 ./ --no-cache
3. make a new container from above image: docker run -it --name test -p 6770:22 svanalysis_starprotocols:1.0.0 /bin/bash
4. login through ssh: ssh root@0.0.0.0 -p 6770
5. you can do anything like in Linux server


Owing to failed to install pacakge 'pysam', we do not install HTSeq and lumpy-sv in our docker image. After you make a container, you can install it in your container using the following commands:

```
#! 1- install miniconda
aria2c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -d ~ && \
bash ~/Miniconda3-latest-Linux-x86_64.sh -b && \
~/miniconda3/bin/conda config --add channels defaults    && \
~/miniconda3/bin/conda config --add channels bioconda    && \
~/miniconda3/bin/conda config --add channels conda-forge && \
rm -f ~/Miniconda3-latest-Linux-x86_64.sh && \
~/miniconda3/bin/conda init

#! 2- install software
~/miniconda3/bin/conda install htseq=2.0.2 && \
~/miniconda3/bin/conda create -n python2.7.15 python=2.7.15 && \
source activate python2.7.15 && \
~/miniconda3/bin/conda install lumpy-sv=0.2.13 && \
conda deactivate 
```

## How to install all software by yourself
```
#! 1- download binary executable files
mkdir -p ~/software && \
aria2c https://github.com/broadinstitute/gatk/archive/refs/tags/4.4.0.0.tar.gz -d ~/software && \
cd ~/software && tar -zxvf gatk-4.4.0.0.tar.gz && rm gatk-4.4.0.0.tar.gz && \
aria2c https://sourceforge.net/projects/svseq2/files/SVseq2_2/SVseq2_2/download -d ~/software/SVseq2 && \
cd ~/software/SVseq2 && chmod 770 ./SVseq2_2 && \
aria2c https://github.com/arq5x/bedtools2/releases/download/v2.30.0/bedtools.static.binary -d ~/software/bedtools_v2.30.0 && \
cd ~/software/bedtools_v2.30.0 && \
mv bedtools.static.binary bedtools && \
chmod a+x bedtools && \
aria2c https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20230116.zip -d ~/software/plink_1.9 && \
cd ~/software/plink_1.9 && unzip plink_linux_x86_64_20230116.zip && \
rm ~/software/plink_1.9/plink_linux_x86_64_20230116.zip && \
aria2c https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20230516.zip -d ~/software/plink_2 && \
cd ~/software/plink_2 && unzip plink2_linux_x86_64_20230516.zip && \
rm ~/software/plink_2/plink2_linux_x86_64_20230516.zip && \
aria2c https://github.com/broadinstitute/picard/releases/download/2.26.11/picard.jar -d ~/software/picard_2.26.11 && \
aria2c https://github.com/DecodeGenetics/svimmer/archive/refs/heads/master.zip -d ~/software && unzip ~/software/svimmer-master.zip && rm ~/software/svimmer-master.zip && \
aria2c https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.4/graphtyper -d ~/software/graphtyper_2.7.4 && \
cd ~/software/graphtyper_2.7.4 && chmod a+x graphtyper

#! 2- install miniconda
aria2c https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -d ~ && \
bash ~/Miniconda3-latest-Linux-x86_64.sh -b && \
~/miniconda3/bin/conda config --add channels defaults    && \
~/miniconda3/bin/conda config --add channels bioconda    && \
~/miniconda3/bin/conda config --add channels conda-forge && \
rm -f ~/Miniconda3-latest-Linux-x86_64.sh && \
~/miniconda3/bin/conda init

#! 3- install other software
~/miniconda3/bin/conda install fastqc=0.12.1 && \
~/miniconda3/bin/conda install bwa=0.7.17 && \
~/miniconda3/bin/conda install samtools=1.17 && \
~/miniconda3/bin/conda install vcftools=0.1.16 && \
~/miniconda3/bin/conda install hisat2=2.2.1 && \
~/miniconda3/bin/conda install trim-galore && \
~/miniconda3/bin/conda install htseq=2.0.2 && \
~/miniconda3/bin/conda install -c bioconda vcflib && \
~/miniconda3/bin/conda install -c bioconda gfold=1.1.4 && \
~/miniconda3/bin/conda install -c bioconda graphtyper=2.7.2
~/miniconda3/bin/conda install repeatmodeler=2.0.1
~/miniconda3/bin/conda install repeatmasker=4.0.9_p2
~/miniconda3/bin/conda create -n python2.7.15 python=2.7.15 && \
source activate python2.7.15 && \
~/miniconda3/bin/conda install manta=1.6.0 && \
~/miniconda3/bin/conda install lumpy-sv=0.2.13 && \
conda deactivate
```

