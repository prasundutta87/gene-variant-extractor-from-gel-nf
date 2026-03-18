FROM rocker/r-base:4.5.2

LABEL maintainer="Prasun Dutta"
LABEL description="Software environment for gene variant extraction pipeline"

RUN apt-get update && apt-get install -y --no-install-recommends \
    # core tools
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
    \
    # bioinformatics tools
    bcftools \
    bedtools \
    samtools \
    tabix \
    \
    # R system dependencies (critical)
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    \
    # font + text rendering (systemfonts / textshaping)
    libfontconfig1-dev \
    libfreetype6-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    \
    # image / graphics stack (ragg)
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libwebp-dev \
    \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install R packages
RUN R -e "options(repos = c(CRAN='https://cloud.r-project.org')); \
          install.packages(c('argparser','duckdb','duckplyr','tidyverse'))"

WORKDIR /workspace

CMD ["/bin/bash"]
