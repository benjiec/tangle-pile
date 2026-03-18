# Pile

Piles of reads, piles of transcripts, piles of alignments, piles of evidence...

This repository contains scripts and pipeline for processing NGS and RNAseq
data. The aim here is to standardize naming and file location conventions for
working with sequences, transcriptomes, alignments across species and samples.


## Setup

Best to use the Docker setup. See Dockerfile for dependencies. See
docker-compose.yml file for mounting local "workspaces" and "reference"
directories onto Docker containers.

There should be a directory where all the workspace and reference files are
stored. Set the following variables on the host, then use `docker-compose` to
run commands.

```
PILE_WORKSPACES_DIR=...
PILE_REFERENCE_DIR=...
```


## Data Model and Directory Setup

Pile works with the following data types

  * reference: files associated with a reference genome, each with a NCBI accession

  * workspace
      * reads: NGS sequencing reads, named with SRA accessions or other unique identifiers
      * transcriptomes: uniquely named transcriptomes, each with .fna, .cds.fna, .faa
      * alignments: sample reads aligned against a transcriptome, named as <sra>-<transcriptome>.bam
      * transcripts: a transcript's alignments from those of a sample against a transcriptome

For more details on the file naming conventions, see the `pile.Defaults` class.


## Docker Commands

All commands require a workspace set using the PILE_WORKSPACE environment
variable, which is passed to the container via `docker-compose.yml`.

Fetching SRA reads

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/sra-fetch.py \
  SRR9331959
```

Indexing a transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/transcriptome-index.py \
  SRR9331959_algae_denovo
```

Align a sample against a transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/transcriptome-align.py \
  SRR9331961 SRR9331959_algae_denovo \
  --cpus 6
```

Extract a transcript's alignments from the alignments of a sample's reads
against a transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/alignment-extract.py -a \
  SRR9331961 SRR9331959_algae_denovo TRINITY_DN7562_c0_g1_i1
``` 


## Insider a Docker Container

The following commands can be run inside a docker container, in case you need to do some custom analysis. To get into a Docker container, run

```
docker-compose run --rm pile /bin/bash
```

Search for a sequence

```
bowtie2 --local -f -p 8 \
  -x /data/doi:10.1126_sciadv.aba2498/transcriptomes/SRR9331959_algae_denovo/transcript_clusters.fna \
  -U query.fasta \
  -S results.sam
```


## Download data from NCBI SRA

Use the following to get a list of SRR accessions to download

```
pip3 install pysradb
pysradb metadata PRJNA591730 > metadata.tsv
```
