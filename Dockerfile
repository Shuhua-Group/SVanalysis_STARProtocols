FROM amd64/ubuntu:focal AS builder

ENV LANG=C.UTF-8 TZ=Asia/Shanghai PS=svprotocol port=22
#设置时区、更新镜像软件、安装aria2（下载工具替代wget，curl以获取更快的下载速度，容错/下载会自动重试）
#openssh服务并更新配置文件，使root账户可以登录、更新root账户密码为设置值
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && \
    echo $TZ > /etc/timezone && \
    apt update && apt upgrade -y && \
    apt install --no-install-recommends openssh-server -y && \
    apt autoremove -y && rm -rf /var/lib/apt/lists/* && \
    sed -i "s/#PermitEmptyPasswords no/PermitEmptyPasswords no/g" /etc/ssh/sshd_config && \
    sed -i "s/#PermitRootLogin prohibit-password/PermitRootLogin yes/g" /etc/ssh/sshd_config && \
    sed -i "s/#ListenAddress 0.0.0.0/ListenAddress 0.0.0.0/g" /etc/ssh/sshd_config && \
    sed -i "s/#LoginGraceTime 2m/LoginGraceTime 2m/g" /etc/ssh/sshd_config && \
    echo root:${PS} | chpasswd && \
    apt-get update && apt-get upgrade && apt-get -y install \
    aria2 unzip \
    gcc make cmake g++ && \
    apt-get -y install --no-install-recommends \
    libbz2-dev \
    libcurl4-gnutls-dev \
    libncurses5-dev \
    libssl-dev \
    liblzma-dev \
    libicu-dev \
    zlib1g-dev \
    pkg-config && \
    apt-get -y install --no-install-recommends \
    parallel apt-utils && \
    apt-get install -y \
    libvcflib-tools libvcflib-dev && \
    apt-get install --no-install-recommends -y \
    default-jre \
    bwa && \
    apt-get install --no-install-recommends -y \
    python2 python3-pip && \
    mkdir -p ~/software && \
    cd ~/software && \
    aria2c https://github.com/DaehwanKimLab/hisat2/archive/refs/heads/master.zip -d ~/software && \
    unzip hisat2-master.zip && \
    rm hisat2-master.zip && \
    cd hisat2-master && \
    make && \
    cd ~/software && \
    aria2c https://github.com/samtools/samtools/releases/download/1.17/samtools-1.17.tar.bz2 -d ~/software && \
    tar -xjf samtools-1.17.tar.bz2 && \
    rm samtools-1.17.tar.bz2 && \
    cd samtools-1.17 && \
    ./configure && \
    make && \
    make install && \
    cp samtools /usr/local/bin && \
    cd htslib-1.17 && make && \
    cd ~/software && \
    aria2c https://github.com/vcftools/vcftools/releases/download/v0.1.16/vcftools-0.1.16.tar.gz -d ~/software && \
    tar -zxvf vcftools-0.1.16.tar.gz && \
    rm vcftools-0.1.16.tar.gz && \
    cd vcftools-0.1.16 && \
    ./configure && \
    make && \
    make install && \
    cd .. && \
    rm -r vcftools-0.1.16 && \
    cd ~/software && aria2c https://github.com/samtools/bcftools/releases/download/1.17/bcftools-1.17.tar.bz2 -d ~/software && \
    tar -vxjf bcftools-1.17.tar.bz2 && \
    rm bcftools-1.17.tar.bz2 && \
    cd bcftools-1.17 && \
    make && \
    make install && \
    cp bcftools /usr/local/bin && \
    cd ../ && rm -rf bcftools-1.17 && \
    cd ~/software && aria2c https://github.com/biod/sambamba/releases/download/v0.8.2/sambamba-0.8.2-linux-amd64-static.gz -d ~/software && \
    gzip -d sambamba-0.8.2-linux-amd64-static.gz && \
    chmod a+x sambamba-0.8.2-linux-amd64-static && \
    mv sambamba /usr/local/bin/. && \
    cd ~/software && \
    aria2c https://github.com/GregoryFaust/samblaster/releases/download/v.0.1.26/samblaster-v.0.1.26.tar.gz -d ~/software && \
    tar -zxvf samblaster-v.0.1.26.tar.gz && \
    rm samblaster-v.0.1.26.tar.gz && \
    cd samblaster-v.0.1.26 && \
    make && \
    mv samblaster /usr/local/bin/. && \
    cd ../ && rm -rf samblaster-v.0.1.26 && \
    ln -s /usr/bin/python2 /usr/bin/python && \
    cd ~/software && aria2c https://github.com/arq5x/lumpy-sv/archive/refs/tags/0.2.13.zip -d ~/software && \
    unzip lumpy-sv-0.2.13.zip && \
    rm lumpy-sv-0.2.13.zip && \
    cd lumpy-sv-0.2.13 && \
    make && \
    cd ~/software && \
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
    cd ~/software && \
    aria2c https://github.com/DecodeGenetics/svimmer/archive/refs/heads/master.zip -d ~/software && unzip ~/software/svimmer-master.zip && rm ~/software/svimmer-master.zip && \
    aria2c https://github.com/DecodeGenetics/graphtyper/releases/download/v2.7.4/graphtyper -d ~/software/graphtyper_2.7.4 && \
    cd ~/software/graphtyper_2.7.4 && chmod a+x graphtyper && \
    aria2c https://github.com/Illumina/manta/releases/download/v1.6.0/manta-1.6.0.centos6_x86_64.tar.bz2 -d ~/software && \
    cd ~/software && \
    tar -xjf manta-1.6.0.centos6_x86_64.tar.bz2 && \
    rm manta-1.6.0.centos6_x86_64.tar.bz2 && \
    # apt-get install --no-install-recommends gnupg -y && \
    # apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && \
    # echo "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/" >> /etc/apt/sources.list && \
    apt-get update && apt-get install r-base -y && \
    apt-get install --no-install-recommends libmariadbclient-dev libpq-dev -y && \
    # apt-get install --no-install-recommends build-essential libxml2-dev libfontconfig1-dev -y && \
    # apt-get install --no-install-recommends -y \
    # libharfbuzz-dev libfribidi-dev  libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev gfortran libmariadbclient-dev libpq-dev && \
    Rscript -e 'install.packages("BiocManager")' && \
    Rscript -e 'install.packages("https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz",repos = NULL,type = "source")' && \
    Rscript -e 'BiocManager::install("edgeR")' && \
    Rscript -e 'install.packages("dplyr",dependencies=TRUE)' && \
    apt remove gcc make cmake g++ -y && \
    apt autoremove -y && rm -rf /var/lib/apt/lists/*
# the above ok 
#ssh port
EXPOSE      $port
# initial image, set password by user, and start ssh service
ENTRYPOINT  ["/bin/bash", "-c", "echo root:${PS} | chpasswd && service ssh start -D"]
