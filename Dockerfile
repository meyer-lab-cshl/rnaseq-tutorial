FROM condaforge/miniforge3:latest

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        less nano vim file tree curl unzip procps htop \
    && rm -rf /var/lib/apt/lists/*

RUN mamba install -y -n base -c conda-forge -c bioconda -c nodefaults \
        snakemake=9.27.1 graphviz \
    && mamba clean -afy

# Snakemake keys each env by its file content and the prefix, so envs built
# here are reused at runtime as long as envs/*.yaml stay unchanged.
ENV SNAKEMAKE_PROFILE=/etc/snakemake
RUN mkdir -p /etc/snakemake \
    && printf 'use-conda: true\nconda-prefix: /opt/snakemake-envs\n' \
        > /etc/snakemake/config.yaml

WORKDIR /build
COPY setup.smk ./
COPY envs ./envs
RUN printf 'rule gffutils:\n    conda: "envs/gffutils.yaml"\n' > extra.smk \
    && printf '    output: "x"\n    shell: "true"\n' >> extra.smk \
    && snakemake -s setup.smk --cores 1 --conda-create-envs-only \
    && snakemake -s extra.smk --cores 1 --conda-create-envs-only \
    && mamba clean -afy && rm -rf /build

WORKDIR /work
CMD ["bash"]
