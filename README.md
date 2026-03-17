# Pile

Piles of reads, piles of transcripts, piles of alignments, piles of evidence...

This repository contains scripts and pipeline for processing NGS and RNAseq
data. The aim here is to standardize naming and file location conventions for
working with sequences, transcriptomes, alignments across species and samples.


## Setup

Best to use the Docker setup. See Dockerfile for dependencies. See
docker-compose.yml file for mounting local "workspaces" and "reference"
directories onto Docker containers.


## Data Model and Directory Setup

Pile works with the following data types

  * reference: files associated with a reference genome, each with a NCBI accession

  * workspace
      * reads: NGS sequencing reads, named with SRA accessions or other unique identifiers
      * transcriptomes: uniquely named transcriptomes, each with .fna, .cds.fna, .faa
      * alignments: sample reads aligned against a transcriptome, named as <sra>-<transcriptome>.{bam}

For more details on the file naming conventions, see `pile.Defaults` class.
