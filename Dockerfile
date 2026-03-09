FROM rocker/r-base:4.5.2

LABEL maintainer="Prasun Dutta"
LABEL description="Software environment for gene variant extraction pipeline"

# Install bioinformatics and bash tools (latest available from Debian repos) and do not install optional suggested packages automatically
#Then finally remove downloaded .deb files and apt metadata to keep the image small
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
    bcftools \
    bedtools \
    samtools \
    tabix \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# Install latest CRAN versions of R packages
RUN R -e "install.packages(c('tidyverse','argparser','duckdb','duckplyr'), repos='https://cloud.r-project.org')"

# Set the default directory
WORKDIR /workspace

# Default command: start an interactive Bash shell when the container runs
# This allows the container to behave like a standard Linux environment
# and makes it easy to execute bash scripts or use it inside Nextflow
CMD ["/bin/bash"]
