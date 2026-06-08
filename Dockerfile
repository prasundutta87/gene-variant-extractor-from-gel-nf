FROM rocker/r-base:4.5.2
LABEL maintainer="Prasun Dutta"
LABEL description="Software environment for gene variant extraction pipeline"

RUN apt-get update && apt-get install -y --no-install-recommends \
    bash \
    coreutils \
    findutils \
    grep \
    sed \
    gawk \
    tar \
    gzip \
    bzip2 \
    less \
    vim \
    wget \
    curl \
    git \
    procps \
    ca-certificates \
    build-essential \
    pkg-config \
    autoconf \
    automake \
    make \
    gcc \
    perl \
    \
    # bioinformatics tools (NOT bcftools — compiled from source below)
    bedtools \
    samtools \
    tabix \
    \
    # htslib/bcftools build dependencies (from official bcftools INSTALL docs)
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libperl-dev \
    libgsl0-dev \
    \
    # R system dependencies
    libxml2-dev \
    \
    # font + text rendering
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    \
    # image / graphics stack
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev \
    \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Build htslib first with explicit libcurl support
ARG HTSLIB_VERSION=1.21
ARG BCFTOOLS_VERSION=1.21

RUN wget https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2

# Build bcftools against the htslib we just installed
RUN wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix=/usr/local --with-htslib=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2

# Print full version output for debugging — does not fail build
RUN bcftools --version

# Install R packages
RUN R -e "options(repos = c(CRAN='https://cloud.r-project.org')); \
          install.packages(c('argparser','duckdb','duckplyr','tidyverse','data.table','gt','ggpubr'))"

WORKDIR /workspace
CMD ["/bin/bash"]