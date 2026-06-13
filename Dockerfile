FROM rocker/tidyverse:4.5.2
# Base image rocker/tidyverse already includes:
#   R, tidyverse, ggplot2, dplyr, gcc, make, perl, git, wget,
#   ca-certificates, zlib1g-dev, libbz2-dev, liblzma-dev, libssl-dev,
#   libcurl4-openssl-dev, libxml2-dev, font and graphics libraries and more.

LABEL maintainer="Prasun Dutta"
LABEL description="Software environment for gene variant extraction pipeline"

# ── System tools and bioinformatics dependencies ─────────────────────────────
# bzip2        : needed to extract .tar.bz2 source archives (not in base image)
# autoconf/automake : needed to run ./configure for htslib/bcftools compilation
# bedtools     : genome arithmetic (intersect, merge, etc.)
# samtools     : SAM/BAM manipulation
# tabix        : index and query genomic files
# libperl-dev  : Perl bindings for bcftools (enables -i/-e filtering with Perl)
# libgsl0-dev  : GNU Scientific Library for bcftools polysomy command
# Note: libcurl4-openssl-dev is already in rocker/tidyverse and provides
#       the libcurl support needed for bcftools S3 streaming
RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    autoconf \
    automake \
    less \
    vim \
    bedtools \
    samtools \
    tabix \
    libperl-dev \
    libgsl0-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# ── bcftools compiled from source with S3/libcurl support ────────────────────
# The apt version of bcftools does NOT include S3 streaming support.
# We compile both htslib and bcftools from source so bcftools can read
# directly from s3:// paths without downloading full VCF shards.
# Steps:
#   1. Download and compile htslib — ./configure auto-detects libcurl and
#      enables S3 support
#   2. Download and compile bcftools against that htslib
#   3. Clean up source files to keep the image small
ARG HTSLIB_VERSION=1.21
ARG BCFTOOLS_VERSION=1.21

RUN wget -q https://github.com/samtools/htslib/releases/download/${HTSLIB_VERSION}/htslib-${HTSLIB_VERSION}.tar.bz2 \
    && tar -xjf htslib-${HTSLIB_VERSION}.tar.bz2 \
    && cd htslib-${HTSLIB_VERSION} \
    && ./configure --prefix=/usr/local \
    && make -j$(nproc) \
    && make install \
    && ldconfig \
    && cd .. \
    && rm -rf htslib-${HTSLIB_VERSION} htslib-${HTSLIB_VERSION}.tar.bz2 \
    && wget -q https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2 \
    && cd bcftools-${BCFTOOLS_VERSION} \
    && ./configure --prefix=/usr/local --with-htslib=/usr/local \
    && make -j$(nproc) \
    && make install \
    && cd .. \
    && rm -rf bcftools-${BCFTOOLS_VERSION} bcftools-${BCFTOOLS_VERSION}.tar.bz2

# ── Additional R packages via Posit Public Package Manager ───────────────────
# Posit PPM provides pre-compiled binary packages for Ubuntu (noble = 24.04).
# This is much faster than compiling from source. Same packages and versions 
# as CRAN, just pre-built for Linux.
RUN R -e "options(repos = c(CRAN='https://packagemanager.posit.co/cran/__linux__/noble/latest')); \
          install.packages(c('argparser','duckdb','duckplyr','data.table','gt','ggpubr','arrow'))"

# ── Verify all tools ──────────────────────────────────────────────────────────
# Build fails here if anything is missing — better to know now than later.
RUN echo "=== Checking all tools ===" && \
    bcftools --version | head -1 && \
    samtools --version | head -1 && \
    bedtools --version && \
    tabix --version 2>&1 | head -1 && \
    Rscript -e "library(tidyverse); library(duckplyr); library(data.table); \
                library(argparser); library(gt); library(ggpubr); library(arrow) \
                cat('All R packages OK\n')" && \
    echo "=== All tools verified ==="

WORKDIR /workspace
CMD ["/bin/bash"]
