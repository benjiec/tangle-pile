import os
import tempfile
import unittest

from pile.transdecoder_to_detected import parse_transdecoder_gff, results_to_detected_table


class TestGFFParsing(unittest.TestCase):

    def test_parses_normal_transdecoder_gff(self):

        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	5443	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1

tx2	transdecoder	gene	1	5491	.	+	.	ID=g2;Name="ORF type:5prime_partial (-)"
tx2	transdecoder	mRNA	1	5491	.	+	.	ID=m2;Parent=g2;Name="ORF type:5prime_partial (-)"
tx2	transdecoder	exon	1	5491	.	+	.	ID=m2.exon;Parent=m2
tx2	transdecoder	CDS	251	5491	.	+	0	ID=cds.m2;Parent=m2;5_prime_partial=true
tx2	transdecoder	five_prime_UTR	1	250	.	+	.	ID=m2.utr5p1;Parent=m2
"""
                f.write(gff)

            results = parse_transdecoder_gff(gff_fn)
            self.assertEqual(len(results), 2)

            self.assertEqual(results[0]["transcript_accession"], "tx1")
            self.assertEqual(results[0]["mrna_accession"], "m1")
            self.assertEqual(results[0]["protein_length"], 1731)
            # reverse strand coordinates
            self.assertEqual(results[0]["strand"], "-")
            self.assertEqual(results[0]["transcript_gene_start"], 5443)
            self.assertEqual(results[0]["transcript_gene_end"], 1)
            self.assertEqual(results[0]["transcript_protein_start"], 5443)
            self.assertEqual(results[0]["transcript_protein_end"], 251)

            self.assertEqual(results[1]["transcript_accession"], "tx2")
            self.assertEqual(results[1]["mrna_accession"], "m2")
            self.assertEqual(results[1]["protein_length"], 1747)
            # reverse strand coordinates
            self.assertEqual(results[1]["strand"], "+")
            self.assertEqual(results[1]["transcript_gene_start"], 1)
            self.assertEqual(results[1]["transcript_gene_end"], 5491)
            self.assertEqual(results[1]["transcript_protein_start"], 251)
            self.assertEqual(results[1]["transcript_protein_end"], 5491)

    def test_flag_multiple_exons(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	exon	3000	5443	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Cannot parse multiple exon entries")

    def test_flag_multiple_cds(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	100	200	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Cannot parse multiple CDS entries")

    def test_flag_mismatched_mrna(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx3	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Mismatched transcript accession")

    def test_flag_mismatched_cds(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx3	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Mismatched transcript accession")

    def test_flag_mismatched_mrna_coordinates(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	2	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Inconsistent mRNA coordinates")

    def test_flag_inconsistent_cds_strand(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	+	0	ID=cds.m1;Parent=m1;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Inconsistent CDS strand")

    def test_flag_inconsistent_cds_parent(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            gff_fn = os.path.join(temp_dir, "proteins.gff")
            with open(gff_fn, 'w') as f:

                gff = """
tx1	transdecoder	gene	1	5443	.	-	.	ID=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	mRNA	1	5443	.	-	.	ID=m1;Parent=g1;Name="ORF type:5prime_partial (-)"
tx1	transdecoder	exon	1	2000	.	-	.	ID=m1.exon;Parent=m1
tx1	transdecoder	CDS	251	5443	.	-	0	ID=cds.m1;Parent=m2;5_prime_partial=true
tx1	transdecoder	three_prime_UTR	1	250	.	-	.	ID=m1.utr3p1;Parent=m1
"""
                f.write(gff)

            with self.assertRaises(AssertionError) as e:
                parse_transdecoder_gff(gff_fn)
            self.assertEqual(str(e.exception), "Mismatched CDS parent accession")


class TestDetectedTSV(unittest.TestCase):

    def test_produces_correct_detection_file_with_protein_coordinates(self):
        results = [
            dict(transcript_accession="tx1",
                 transcript_gene_start=5443,     # reverse strand gene
                 transcript_gene_end=1,
                 gene_accession="g1",
                 mrna_accession="m1",
                 transcript_protein_start=5440,  # reverse strand protein - will need to re-orient with gene direction
                 transcript_protein_end=251,
                 protein_length=1730,
                 strand="-"),
            dict(transcript_accession="tx2",
                 transcript_gene_start=1,        # forward strand gene
                 transcript_gene_end=5491,
                 gene_accession="g2",
                 mrna_accession="m2",
                 transcript_protein_start=251,   # forward strand protein
                 transcript_protein_end=5488,
                 protein_length=1746,
                 strand="+")
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            detected_tsv = os.path.join(temp_dir, "detected.tsv")
            results_to_detected_table(results, detected_tsv, "foo", batch="20260327_ff30c5d5")

            expected = """
detection_type	detection_method	batch	query_accession	query_database	query_type	target_accession	target_database	target_type	target_model	query_start	query_end	target_start	target_end	evalue	bitscore	bitscore_threshold	custom_metric_name	custom_metric_value
sequence	transdecoder	20260327_ff30c5d5	tx1	foo	transcript	g1	foo	gene		5443	1	1	5443					
sequence	transdecoder	20260327_ff30c5d5	g1	foo	gene	m1	foo	protein		4	5193	1	1730					
sequence	transdecoder	20260327_ff30c5d5	tx2	foo	transcript	g2	foo	gene		1	5491	1	5491					
sequence	transdecoder	20260327_ff30c5d5	g2	foo	gene	m2	foo	protein		251	5488	1	1746					
"""

            with open(detected_tsv, "r") as f:
                self.maxDiff = None
                self.assertEqual(f.read().strip(), expected.strip())

    def test_flag_non_unique_gene(self):
        results = [
            dict(transcript_accession="tx1",
                 transcript_gene_start=5443,
                 transcript_gene_end=1,
                 gene_accession="g1",
                 mrna_accession="m1",
                 transcript_protein_start=5443,
                 transcript_protein_end=251,
                 protein_length=1731,
                 strand="-"),
            dict(transcript_accession="tx2",
                 transcript_gene_start=1,
                 transcript_gene_end=5491,
                 gene_accession="g1",
                 mrna_accession="m2",
                 transcript_protein_start=251,
                 transcript_protein_end=5491,
                 protein_length=1747,
                 strand="+")
        ]

        with tempfile.TemporaryDirectory() as temp_dir:
            detected_tsv = os.path.join(temp_dir, "detected.tsv")
            with self.assertRaises(AssertionError) as e:
                results_to_detected_table(results, detected_tsv, "foo", batch="20260327_ff30c5d5")
            self.assertEqual(str(e.exception), "Duplicate gene ID g1")
