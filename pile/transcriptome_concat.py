#!/usr/bin/env python3

from pile import run_command
from pile.defaults import Defaults
from tangle.sequence import read_fasta_as_dict, write_fasta_from_dict


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptomes", nargs='+')
    parser.add_argument("--unclustered", action="store_true", default=False)
    args = parser.parse_args()

    fns = []
    for tr in args.transcriptomes:
        if args.unclustered:
            fn = Defaults.transcriptome_unclustered_fasta(Defaults.workspace(), tr)
        else:
            fn = Defaults.transcriptome_fasta(Defaults.workspace(), tr)
        fns.append(fn)

    transcriptome_name = "_".join(args.transcriptomes)
    if args.unclustered:
        target_fn = Defaults.transcriptome_unclustered_fasta(Defaults.workspace(), transcriptome_name, gzip=True)
    else:
        target_fn = Defaults.transcriptome_fasta(Defaults.workspace(), transcriptome_name, gzip=True)

    sequences = {}
    for fn in fns:
        print(fn)
        d = read_fasta_as_dict(fn)

        old_count = len(sequences)
        add_count = len(d)
        sequences.update(d)
        new_count = len(sequences)
        if new_count != old_count+add_count:
            print(f"  new count of sequences = {new_count}, different from expected {old_count+add_count}, sequence accessions may not be unique")

    print(target_fn)
    write_fasta_from_dict(sequences, target_fn)
