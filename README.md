# Pile

Piles of reads, piles of transcripts, piles of alignments, piles of evidence...

This repository contains scripts and pipeline for processing NGS and RNAseq
data.


## Setup

MMSeqs2 docker image: `docker pull ghcr.io/soedinglab/mmseqs2`

Conda environment

```
conda create -n pile_env -c bioconda -c conda-forge
conda activate pile_env
conda install -c bioconda bbmap salmon bowtie2 samtools transdecoder trinity
```

Install NCBI SRA tools, see https://github.com/ncbi/sra-tools


## Data Model and Directory Setup

Pile works with the following data types

  * reference: files associated with a reference genome, each with a NCBI accession

  * workspace
      * reads: NGS sequencing reads, named with SRA accessions or other unique identifiers
      * transcriptome: uniquely named transcriptomes, e.g. .fna, .cds.fna, .faa
      * sample-alignments: a sample's reads, aligned against a transcriptome
