FROM mambaorg/micromamba:latest

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

ARG MAMBA_DOCKERFILE_ACTIVATE=1

# 2. Set the working directory inside the container
WORKDIR /app

# 3. Copy your requirements and install them
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copy your actual script/logic into the image
COPY pile ./pile

# (Optional) Add your app directory to the PATH
ENV PATH="/app:${PATH}"
ENV PYTHONPATH="/app"

# We don't use ENTRYPOINT or CMD because Nextflow 
# will override them to run its own wrapper script.
