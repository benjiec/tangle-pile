from tangle import open_file_to_read, unique_batch
from tangle.detected import DetectedTable
from pile import Path
from pile.defaults import Defaults


def parse_transdecoder_gff(file_path):
    results = []
    current_entry = {}

    with open_file_to_read(file_path) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            cols = line.strip().split('\t')
            feature_type = cols[2]
            attributes = dict(item.split('=') for item in cols[8].split(';') if '=' in item)
            
            # State: Start of a new gene block
            if feature_type == 'gene':
                current_entry = {
                    'transcript_accession': cols[0],
                    'gene_accession': attributes.get('ID'),
                    'strand': cols[6]
                }
                if current_entry["strand"] == "+":
                    current_entry.update({
                        'transcript_gene_start': int(cols[3]),
                        'transcript_gene_end': int(cols[4])
                    })
                elif current_entry["strand"] == "-":
                    current_entry.update({
                        'transcript_gene_start': int(cols[4]),
                        'transcript_gene_end': int(cols[3])
                    })

            elif feature_type == 'mRNA':
                assert cols[0] == current_entry["transcript_accession"], "Mismatched transcript accession"
                assert sorted([int(cols[3]), int(cols[4])]) == sorted([current_entry['transcript_gene_start'], current_entry['transcript_gene_end']]), "Inconsistent mRNA coordinates"
                current_entry['mrna_accession'] = attributes.get('ID')
                current_entry["exon_count"] = 0

            elif feature_type == 'exon':
                assert cols[0] == current_entry["transcript_accession"], "Mismatched transcript accession"
                assert current_entry["exon_count"] == 0, "Cannot parse multiple exon entries"
                current_entry["exon_count"] += 1

            elif feature_type == 'CDS':
                raw_cds_id = attributes.get('ID', '')
                parent_id = attributes.get('Parent', '')
                
                assert cols[0] == current_entry["transcript_accession"], "Mismatched transcript accession"
                assert parent_id == current_entry['mrna_accession'], "Mismatched CDS parent accession"
                assert cols[6] == current_entry['strand'], "Inconsistent CDS strand"
                assert "protein_length" not in current_entry, "Cannot parse multiple CDS entries"

                if current_entry["strand"] == "+":
                    current_entry.update({
                        'transcript_protein_start': int(cols[3]),
                        'transcript_protein_end': int(cols[4])
                    })
                elif current_entry["strand"] == "-":
                    current_entry.update({
                        'transcript_protein_start': int(cols[4]),
                        'transcript_protein_end': int(cols[3])
                    })

                current_entry["protein_length"] = (abs(int(cols[3]) - int(cols[4])) + 1) // 3
                results.append(current_entry)
                
    return results


def results_to_detected_table(
    results,
    result_tsv,
    transcriptome_name,
    batch=None
  ):

    rows = []
    batch = unique_batch() if batch is None else batch
    unique_gene_ids = {}

    for res in results:
        assert res["gene_accession"] not in unique_gene_ids, f"Duplicate gene ID {res['gene_accession']}"
        unique_gene_ids[res["gene_accession"]] = 1

    for res in results:
        row = {}
        row["detection_type"] = "sequence"
        row["detection_method"] = "transdecoder"
        row["batch"] = batch
        row["query_database"] = transcriptome_name
        row["target_database"] = transcriptome_name

        gene_row = row.copy()
        gene_row["query_type"] = "transcript"
        gene_row["target_type"] = "gene"
        gene_row["query_accession"] = res["transcript_accession"]
        gene_row["target_accession"] = res["gene_accession"]
        gene_row["query_start"] = res["transcript_gene_start"]
        gene_row["query_end"] = res["transcript_gene_end"]
        gene_row["target_start"] = 1
        gene_row["target_end"] = abs(res["transcript_gene_start"]-res["transcript_gene_end"])+1

        protein_row = row.copy()
        protein_row["query_type"] = "gene"
        protein_row["target_type"] = "protein"
        protein_row["query_accession"] = res["gene_accession"]
        protein_row["target_accession"] = res["mrna_accession"]
        protein_row["target_start"] = 1
        protein_row["target_end"] = res["protein_length"]

	# protein start coordinates, inside gene, is relative to the gene
	# orientation, not transcript orientation. in gff, the coordinate
	# system is relative to the contig/sequence, and is independent of the
	# hierarchy of features. in the detected table format, the coordinate
	# system is part of the hierarchy of features, so we need to make
	# adjustments.

        if res["strand"] == "+":
            offset = res["transcript_gene_start"] - 1
            protein_row["query_start"] = res["transcript_protein_start"] - offset
            protein_row["query_end"] = res["transcript_protein_end"] - offset

        elif res["strand"] == "-":
            gene_length = abs(res["transcript_gene_start"]-res["transcript_gene_end"])+1
            protein_row["query_start"] = (gene_length-res["transcript_protein_start"])+1
            protein_row["query_end"] = (gene_length-res["transcript_protein_end"])+1

        rows.append(gene_row)
        rows.append(protein_row)

    DetectedTable.write_tsv(result_tsv, rows)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("transcriptome")
    args = parser.parse_args()

    gff_file = Defaults.transcriptome_proteins_gff(Defaults.workspace(), args.transcriptome)
    transcriptome_dir = Defaults.transcriptome_dir(Defaults.workspace(), args.transcriptome)
    detected_tsv = str(Path(transcriptome_dir) / "detected.tsv")

    results = parse_transdecoder_gff(gff_file)
    results_to_detected_table(results, detected_tsv, args.transcriptome)
