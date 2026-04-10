#!/usr/bin/env python3

from pathlib import Path
from pile import run_command
from pile.defaults import Defaults


def bowtie2_filter_reads_remove(read_1, read_2, genomic_fn, prefix, cpus):

    # not using --local because end-to-end alignment is more stringent and we
    # are keeping all reads that don't map fully to the specified genome

    run_command(
      "bowtie2", "-p", str(cpus),
      "-x", genomic_fn,
      "-1", reads_1,
      "-2", reads_2,
      "--no-unal",
      "--un-conc", prefix
    )


def bowtie2_filter_reads_capture(read_1, read_2, genomic_fn, prefix, cpus):

    # using --local because we want all reads that map to specified genome even
    # if locally

    run_command(
      "bowtie2", "-p", str(cpus),
      "--local",
      "-x", genomic_fn,
      "-1", reads_1,
      "-2", reads_2,
      "--no-unal",
      "--al-conc", prefix
    )


def normalize_reads(read_1, read_2, normalized_read_1, normalized_read_2, target_depth):
    run_command(
        "bbnorm.sh",
        "-in1="+read_1,
        "-in2="+read_2,
        "-out1="+normalized_read_1,
        "-out2="+normalized_read_2,
        "target="+str(target_depth)
    )


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptome")
    parser.add_argument("sra_accession")
    parser.add_argument("--remove", default=None)
    parser.add_argument("--capture", default=None)
    parser.add_argument("--cpus", default=4, type=int)
    parser.add_argument("--force", default=False, action="store_true")
    args = parser.parse_args()

    orig_reads_1 = reads_1 = Defaults.read_1(Defaults.workspace(), args.sra_accession)
    orig_reads_2 = reads_2 = Defaults.read_2(Defaults.workspace(), args.sra_accession)

    transcriptome_dir = Defaults.transcriptome_dir(Defaults.workspace(), args.transcriptome)
    to_clean = []

    if args.remove:
        genomic_fn = Defaults.ncbi_genome_genomic(args.remove)
        prefix = str(Path(transcriptome_dir) / (args.sra_accession+f"_{args.remove}_removed"))
        output_reads_1 = Path(prefix+".1.fastq")
        output_reads_2 = Path(prefix+".2.fastq")
        if args.force or not output_reads_1.exists() or output_reads_1.stat().st_size == 0:
            bowtie2_filter_reads_remove(reads_1, reads_2, genomic_fn, prefix+".fastq", args.cpus)
        reads_1 = str(output_reads_1)
        reads_2 = str(output_reads_2)
        to_clean.extend([reads_1, reads_2])

    if args.capture:
        genomic_fn = Defaults.ncbi_genome_genomic(args.capture)
        prefix = str(Path(transcriptome_dir) / (args.sra_accession+f"_{args.capture}_captured"))
        output_reads_1 = Path(prefix+".1.fastq")
        output_reads_2 = Path(prefix+".2.fastq")
        if args.force or not output_reads_1.exists() or output_reads_1.stat().st_size == 0:
            bowtie2_filter_reads_capture(reads_1, reads_2, genomic_fn, prefix+".fastq", args.cpus)
        reads_1 = str(output_reads_1)
        reads_2 = str(output_reads_2)
        to_clean.extend([reads_1, reads_2])

    filtered_reads_1 = Defaults.transcriptome_filtered_read_1(Defaults.workspace(), args.transcriptome, args.sra_accession)
    filtered_reads_2 = Defaults.transcriptome_filtered_read_2(Defaults.workspace(), args.transcriptome, args.sra_accession)
    normalize_reads(reads_1, reads_2, filtered_reads_1, filtered_reads_2, 50)

    for fn in to_clean:
        if fn not in (orig_reads_1, orig_reads_2):
            pfn = Path(fn)
            if pfn.exists():
                print(f"removing {fn}")
                pfn.unlink()
