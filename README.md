# Pile

Piles of reads, piles of transcripts, piles of alignments, piles of evidence...

This repository contains scripts and pipeline for processing NGS and RNAseq
data. The aim here is to standardize naming and file location conventions for
working with sequences, transcriptomes, alignments across species and samples.


## Setup

There should be a directory where all the workspace and reference files are
stored. Set the following variables on the host, then use `docker-compose` to
run commands.

```
PILE_WORKSPACES_DIR=...
```

All other commands require a workspace to be set using the PILE_WORKSPACE
environment variable, which is passed to the container via
`docker-compose.yml`.

The following alias is useful

```
alias pile-py='docker-compose -f pile/docker-compose.yml run --rm pile python3'
alias pile-run='docker-compose -f pile/docker-compose.yml run --rm pile'
```

For clustering sequences, get the MMSeqs2 docker image: `docker pull
ghcr.io/soedinglab/mmseqs2`


## Data Model and Directory Setup

Pile works with the following data types

  * reference: files associated with a reference genome, each with a NCBI accession

  * workspace
      * reads: NGS sequencing reads, named with SRA accessions or other unique identifiers
      * transcriptomes: uniquely named transcriptomes, each with .fna, .cds.fna, .faa
      * alignments: sample reads aligned against a transcriptome, named as <sra>-<transcriptome>.bam
      * transcripts: a transcript's alignments from those of a sample against a transcriptome

For more details on the file naming conventions, see the `pile.Defaults` class.


## Quantification workflow

Fetching SRA reads

```
PILE_WORKSPACE=PM32426508 pile-py pile/sra_fetch.py \
  SRR9331959
```

This may be helpful to download all the SRR numbers in bulk

```
pip3 install pysradb
pysradb metadata PRJNA591730 > metadata.tsv
```

Assuming you have a transcriptome -- i.e. a `transcripts.fna` file in the
transcriptome directory, index the transcriptome like this

(NOTE - you probably want to cluster a transcriptome you downloaded from NCBI
or article supplemental data section, see cluster section below)

```
PILE_WORKSPACE=PM32426508 pile-py pile/transcriptome_index.py \
  SRR9331959_algae_denovo
```

Create Salmon index, then use Salmon to map reads to transcriptomes, e.g.

# XXX fix the following and use docker commands
# XXX always put quants in a standard directory with transcriptome

```
# in docker, in transcriptome directory
salmon index -t pseudodiploria_symb.fna.gz_rep_seq.fna -i pseudodiploria_symb.salmon_index
salmon quant -i pseudodiploria_symb.salmon_index \
  -l A --validateMappings -o symb_quants/SRR6255880 -p 2\
  -1 reads/SRR6255880_1.fastq -2 reads/SRR6255880_2.fastq
```


## Building a transcriptome

A transcriptome should have these files in its transcriptome directory

  * `transcripts.fna` - transcripts, likely clustered
  * `proteins.faa` - proteins associated with the transcripts
  * `proteins.gff3` - GFF mapping transcripts to proteins


### Trinity assembly

To filter reads for Trinity assembly

```
PILE_WORKSPACE=PM32426508 pile-py pile/assembly_filter_reads.py \
    SRR9331961_algae_denovo \
    SRR9331961 \
    --remove GCA_014633955.1 \
    --capture GCA_947184155.2
```

Trinity assembly, after you have filtered the reads

```
# inside docker, in transcriptome directory
Trinity --seqType fq --max_memory 20G --CPU 8 --no_normalize_reads \
        --left SRR9331961_filtered_1.fastq \
        --right SRR9331961_filtered_2.fastq \
        --output trinity_transcript 

TrinityStats.pl transcript.Trinity.fasta
```

### Clustering

If you have mmseqs Docker image, then use the following to cluster transcripts.
Note that clustered transcripts should only be used for certain differential
expression analyses where you don't care about alleles and want to maximize
signal. For most genetics work, unclustered transcripts may be best.

```
# on host, in transcriptome directory
<pile_repo_directory>/scripts/mmseqs-cluster-trinity-transcripts \
  transcripts.fna
mv transcripts.fna transcripts.unclustered.fna
mv transcripts.fna_rep_seq.fna.gz transcript.fna.gz
mv transcripts.fna_cluster.tsv transcript_clusters.tsv
```

Once you have clustered, the downstream workflows -- transdecoder,
quantification, classification etc -- all use clustered transcriptome. This
consistency is important.


### TransDecoder ORF and protein prediction

Use TransDecoder to predict ORFs

```
# inside docker, in transcriptome directory
TransDecoder.LongOrfs -t transcripts.fna -m 100
cp transcripts.fna.transdecoder_dir/longest_orfs.pep proteins.faa
cp transcripts.fna.transdecoder_dir/longest_orfs.gff3 proteins.gff3
rm -rf transcripts.fna.transdecoder_dir
```

Converting transcriptome GFF file, from transdecoder, into detected table (you
need to run TransDecoder first)

```
PILE_WORKSPACE=PM32426508 pile-py pile/transdecoder_to_detected.py \
  SRR9331959_algae_denovo
```


## Reviewing/debugging transcripts

Use the following to make a bowtie index, on a path off the top workspaces
directory (i.e. not specific to a workspace). E.g.

```
pile-py pile/bowtie_index.py ncbi/GCA_014633955.1/genomic.fna
```

Align a sample against a transcriptome

```
PILE_WORKSPACE=PM32426508 pile-py pile/transcriptome_align.py \
  SRR9331961 SRR9331959_algae_denovo \
  --cpus 6
```

Extract a specific transcript's alignments from sample alignments (see above)

```
PILE_WORKSPACE=PM32426508 pile-py pile/alignment_extract.py -a \
  SRR9331961 SRR9331959_algae_denovo TRINITY_DN7562_c0_g1_i1
``` 

Create alignment pileups of transcripts in a second transcriptome, releated to
a transcript from a first transcriptome

```
PILE_WORKSPACE=PM32426508 pile-run bash -c \
  "pile/transcriptome_find.py GCA_947184155.2 CAL1161012.1 | \
   pile/transcriptome_search.py SRR9331959_algae_denovo - | \
   pile/alignment_extract.py SRR9331959 SRR9331959_algae_denovo -"
```

Same as above, but find the protein sequence instead of creating pileups

```
PILE_WORKSPACE=PM32426508 pile-run bash -c \
  "pile/transcriptome_find.py GCA_947184155.2 CAL1161012.1 | \
   pile/transcriptome_search.py SRR9331959_algae_denovo - | \
   pile/transcriptome_find.py -f -p SRR9331959_algae_denovo -"
```

Search for a sequence

```
# inside docker
bowtie2 --local -f -p 8 \
  -x /users/pile/PM32426508/transcriptomes/SRR9331959_algae_denovo/transcripts.fna \
  -U query.fasta \
  -S results.sam
```
