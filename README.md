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

For clustering sequences, get the MMSeqs2 docker image: `docker pull ghcr.io/soedinglab/mmseqs2`


## Data Model and Directory Setup

Pile works with the following data types

  * reference: files associated with a reference genome, each with a NCBI accession

  * workspace
      * reads: NGS sequencing reads, named with SRA accessions or other unique identifiers
      * transcriptomes: uniquely named transcriptomes, each with .fna, .cds.fna, .faa
      * alignments: sample reads aligned against a transcriptome, named as <sra>-<transcriptome>.bam
      * transcripts: a transcript's alignments from those of a sample against a transcriptome

For more details on the file naming conventions, see the `pile.Defaults` class.


## Using Docker Image

Use the following to make a bowtie index, on a path off the top workspaces
directory (i.e. not specific to a workspace). E.g.

```
docker-compose run --rm pile \
  python3 pile/bowtie_index.py ncbi/GCA_014633955.1/genomic.fna
```

All other commands require a workspace to be set using the PILE_WORKSPACE
environment variable, which is passed to the container via
`docker-compose.yml`.

Fetching SRA reads

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/sra_fetch.py \
  SRR9331959
```

Indexing a transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/transcriptome_index.py \
  SRR9331959_algae_denovo
```

Align a sample against a transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/transcriptome_align.py \
  SRR9331961 SRR9331959_algae_denovo \
  --cpus 6
```

Extract a transcript's alignments from the alignments of a sample's reads
against a transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/alignment_extract.py -a \
  SRR9331961 SRR9331959_algae_denovo TRINITY_DN7562_c0_g1_i1
``` 

Create alignment pileups of transcripts in a second transcriptome, releated to a transcript from a first transcriptome

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile bash -c \
  "pile/transcriptome_find.py GCA_947184155.2 CAL1161012.1 | \
   pile/transcriptome_search.py SRR9331959_algae_denovo - | \
   pile/alignment_extract.py SRR9331959 SRR9331959_algae_denovo -"
```

Same as above, but find the protein sequence instead of creating pileups

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile bash -c \
  "pile/transcriptome_find.py GCA_947184155.2 CAL1161012.1 | \
   pile/transcriptome_search.py SRR9331959_algae_denovo - | \
   pile/transcriptome_find.py -f -p SRR9331959_algae_denovo -"
```

To filter reads for Trinity assembly

```
PILE_WORKSPACE=doi:10.1126_sciadv.aba2498 docker-compose run --rm pile \
  python3 pile/assembly_filter_reads.py \
    SRR9331961_algae_denovo \
    SRR9331961 \
    --remove GCA_014633955.1 \
    --capture GCA_947184155.2
```

Actual Trinity assembly, do that inside a Docker container since there are
various flags to use and memory issues to monitor.


## Inside Docker Container

The following commands can be run inside a docker container, in case you need to do some custom analysis. To get into a Docker container, run

```
docker-compose run --rm pile /bin/bash
```

Search for a sequence

```
bowtie2 --local -f -p 8 \
  -x /users/pile/doi:10.1126_sciadv.aba2498/transcriptomes/SRR9331959_algae_denovo/transcripts.fna \
  -U query.fasta \
  -S results.sam
```

Use TransDecoder to predict ORFs

```
TransDecoder.LongOrfs -t transcripts.fna -m 100
cp transcripts.fna.transdecoder_dir/longest_orfs.pep proteins.faa
cp transcripts.fna.transdecoder_dir/longest_orfs.gff3 proteins.gff3
rm -rf transcripts.fna.transdecoder_dir
```

Trinity assembly, after you have filtered the reads. First cd to the workspace
transcriptome directory.

```
Trinity --seqType fq --max_memory 20G --CPU 8 --no_normalize_reads \
        --left SRR9331961_filtered_1.fastq \
        --right SRR9331961_filtered_2.fastq \
        --output trinity_transcript 

TrinityStats.pl transcript.Trinity.fasta
```


# Clustering Sequences

If you have mmseqs Docker image, then use the following to cluster transcripts.
Note that transcript clusters should only be used for quantification, when you
don't want to dilute signal. For most genetics work, unclustered, assembled
transcripts may be best.

```
cd <transcriptome_directory>
<pile_repo_directory>/scripts/mmseqs-cluster-trinity-transcripts \
  transcripts.fna
mv transcripts.fna_rep_seq.fna.gz transcript_clusters.fna.gz
mv transcripts.fna_cluster.tsv transcript_clusters.tsv
```


## Download data from NCBI SRA

Use the following to get a list of SRR accessions to download

```
pip3 install pysradb
pysradb metadata PRJNA591730 > metadata.tsv
```
