#!/usr/bin/env python3

import os
import subprocess
from pile import Defaults, process_file_or_literal, FastaHeaderProvenance
from pile.fasta import read_fasta_as_dict


def bowtie2_search_sequence(fn, query_sequences):

    query_text = []
    for a,s in query_sequences.items():
        query_text.append(f">{a}\n{s}\n")
    query_text = "".join(query_text)

    cmd = f"bowtie2 --local -x {fn} -f -U - | samtools view -F 4 | cut -f 1,3,6 | sort -u"

    result = subprocess.run(
        cmd,
        input=query_text,
        shell=True,
        capture_output=True,
        text=True
    )

    results = []
    for l in result.stdout.splitlines():
        results.append(l.split("\t"))
    return results


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptome")
    parser.add_argument("sequence_file")
    parser.add_argument("--no-fasta", action="store_true", default=False)
    args = parser.parse_args()

    tx_fn = Defaults.transcriptome_fasta(Defaults.workspace(), args.transcriptome)
    if not os.path.exists(tx_fn+".1.bt2"):
        raise Exception("Cannot find bowtie index")

    query_sequences = {}
    if args.no_fasta:
        process_file_or_literal(False, args.sequence_file, lambda s: query_sequences.update(f"query_{len(query_sequences.keys())}", s))
    else:
        query_sequences = read_fasta_as_dict(args.sequence_file, preserve_full_accession=True)

    results = bowtie2_search_sequence(tx_fn, query_sequences)
    for query, hit, cigar in results:
        provenance = FastaHeaderProvenance.unpack(query)
        with_cigar = FastaHeaderProvenance.unpack(FastaHeaderProvenance.add_and_str(cigar, provenance))
        print(FastaHeaderProvenance.add_and_str(hit, with_cigar))
