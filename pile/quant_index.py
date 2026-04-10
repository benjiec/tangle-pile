#!/usr/bin/env python3

from pile import run_command, Defaults

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptome")
    parser.add_argument("--unclustered", action="store_true", default=False)
    parser.add_argument("--cpus", default=2)
    args = parser.parse_args()

    if args.unclustered:
        fasta_fn = Defaults.transcriptome_unclustered_fasta(Defaults.workspace(), args.transcriptome)
        index_fn = Defaults.transcriptome_unclustered_salmon_index(Defaults.workspace(), args.transcriptome)
    else:
        fasta_fn = Defaults.transcriptome_fasta(Defaults.workspace(), args.transcriptome)
        index_fn = Defaults.transcriptome_salmon_index(Defaults.workspace(), args.transcriptome)

    run_command(
      "salmon", "index",
      "-t", fasta_fn,
      "-i", index_fn)
