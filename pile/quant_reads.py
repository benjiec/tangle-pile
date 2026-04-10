#!/usr/bin/env python3

from pile import run_command, Defaults

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptome")
    parser.add_argument("sra_accession")
    parser.add_argument("--unclustered", action="store_true", default=False)
    parser.add_argument("--cpus", default=2)
    args = parser.parse_args()

    if args.unclustered:
        index_fn = Defaults.transcriptome_unclustered_salmon_index(Defaults.workspace(), args.transcriptome)
    else:
        index_fn = Defaults.transcriptome_salmon_index(Defaults.workspace(), args.transcriptome)

    quant_dir = Defaults.quant_dir(Defaults.workspace(), args.transcriptome, args.sra_accession)
    read_1 = Defaults.read_1(Defaults.workspace(), args.sra_accession)
    read_2 = Defaults.read_2(Defaults.workspace(), args.sra_accession)

    run_command(
      "salmon", "quant",
      "-l", "A", "--validateMappings", "-p", "2",
      "-i", index_fn,
      "-o", quant_dir,
      "-1", read_1, "-2", read_2)
