#!/usr/bin/env python3

import sys
from pile import Defaults, process_file_or_literal, FastaHeaderProvenance
from pile.fasta import read_fasta_as_dict


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptome")
    parser.add_argument("transcript")
    parser.add_argument("-p", "--protein", action="store_true", default=False, help="Look for a protein accession, by default look for accession in transcript file")
    parser.add_argument("-f", "--file", action="store_true", default=False)
    parser.add_argument("--exact", action="store_true", default=False)
    args = parser.parse_args()

    if args.protein:
        fasta_fn = Defaults.transcriptome_proteins_fasta(Defaults.workspace(), args.transcriptome)
    else:
        fasta_fn = Defaults.transcriptome_fasta(Defaults.workspace(), args.transcriptome)
    entries = read_fasta_as_dict(fasta_fn)

    accessions = []
    process_file_or_literal(not args.file, args.transcript, lambda v: accessions.append(v))

    for fasta_header in accessions:
        provenance = FastaHeaderProvenance.unpack(fasta_header)
        accession = provenance[0]

        if accession in entries:
            print(f">{fasta_header}\n{entries[accession]}")
        elif args.exact:
            raise Exception(f"Cannot find sequence with accession {accession}")
        else:
            found = False
            for k,v in entries.items():
                if "_"+accession.lower()+"_" in k.lower():
                    new_header = FastaHeaderProvenance.add_and_str(k, provenance)
                    print(f">{new_header}\n{v}")
                    found = True
            if not found:
                for k,v in entries.items():
                    if accession.lower() in k.lower() or k.lower() in accession.lower():
                        new_header = FastaHeaderProvenance.add_and_str(k, provenance)
                        print(f">{new_header}\n{v}")
                        found = True
            if not found:
                raise Exception(f"Cannot find sequence with accession {accession}")
