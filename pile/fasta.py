import sys
import gzip
from contextlib import contextmanager


@contextmanager
def open_fasta_to_read(fn):
    if fn == "-":
        yield sys.stdin
    else:
        with open(fn, 'rb') as test_f:
            is_gz = test_f.read(2) == b'\x1f\x8b'
        opener = gzip.open if is_gz else open
        with opener(fn, "rt", encoding="utf-8") as f:
            yield f


def read_fasta_as_dict(path, preserve_full_accession=False):
    sequences_by_accession = {}
    current_acc = None
    current_seq_parts = []

    with open_fasta_to_read(path) as f:
        for raw_line in f:
            if not raw_line:
                continue
            line = raw_line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                # Flush previous
                if current_acc is not None:
                    sequences_by_accession[current_acc] = "".join(current_seq_parts)
                header_content = line[1:].strip()
                if preserve_full_accession:
                    accession = header_content
                else:
                    # accession is the first whitespace-delimited token
                    accession = header_content.split(None, 1)[0]
                current_acc = accession
                current_seq_parts = []
            else:
                current_seq_parts.append(line.strip())
        # Flush final
        if current_acc is not None:
            sequences_by_accession[current_acc] = "".join(current_seq_parts)

    return sequences_by_accession
