#Code for creating the dockerfile

FROM continuumio/miniconda3

USER root

# Pacchetti di sistema generici
RUN apt-get update && apt-get install -y \
    wget unzip gzip curl git && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Installa mamba per gestire gli ambienti conda più rapidamente
RUN conda install -n base -c conda-forge mamba -y

# Copia gli environment YAML forniti
COPY envs/ /envs/

# Crea gli ambienti esistenti
RUN mamba env create -f /envs/sra-env.yml -p /opt/conda/envs/sra && \
    mamba env create -f /envs/qc-env.yml -p /opt/conda/envs/qc && \
    mamba env create -f /envs/align-env.yml -p /opt/conda/envs/align

RUN mamba create -y -p /opt/conda/envs/rnaseq-env -c conda-forge -c bioconda \
    r-base=4.2 \
    bioconductor-deseq2 \
    bioconductor-edger \
    bioconductor-limma \
    bioconductor-tximport \
    bioconductor-apeglm \
    bioconductor-biomart \
    r-pheatmap \
    bioconductor-ihw \
    bioconductor-genomicfeatures \
    bioconductor-annotationdbi \
    bioconductor-org.hs.eg.db \
    bioconductor-org.mm.eg.db \
    bioconductor-clusterprofiler \
    bioconductor-biocparallel \
    bioconductor-complexheatmap \
    r-tidyverse \
    r-data.table \
    r-rmarkdown \
    r-devtools

# Inizializza conda automaticamente in ogni shell
RUN echo "source /opt/conda/etc/profile.d/conda.sh" >> /etc/bash.bashrc

ENTRYPOINT ["bash", "--rcfile", "/etc/bash.bashrc", "-c"]
SHELL ["/bin/bash", "-c"]

# Default: entra in bash con conda attivabile
CMD ["bash", "-c", "source /opt/conda/etc/profile.d/conda.sh && exec bash"]

