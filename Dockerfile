FROM mambaorg/micromamba:latest

USER root
RUN apt-get update && apt-get install -y \
    git

USER $MAMBA_USER
RUN micromamba install -y -n base -c conda-forge -c bioconda \
    python=3.11 \
    pip=26.0.1 \
    bbmap=39.52 \
    salmon=1.10.3 \
    bowtie2=2.5.5 \
    samtools=1.22.1 \
    transdecoder=5.7.1 \
    trinity=2.15.2 \
    sra-tools=3.2.1 \
    && micromamba clean --all --yes

WORKDIR /pile

COPY --chown=$MAMBA_USER:$MAMBA_USER pyproject.toml ./ 
RUN micromamba run -n base pip install --no-cache-dir . || true

COPY --chown=$MAMBA_USER:$MAMBA_USER . .
RUN micromamba run -n base pip install --no-cache-dir .

ENV PATH="/pile:${PATH}"
ENV PYTHONPATH="/pile"

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# We don't use ENTRYPOINT or CMD because Nextflow 
# will override them to run its own wrapper script.
