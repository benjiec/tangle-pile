from pile import run_command, Defaults


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("workspace")
    parser.add_argument("sra_accession")
    parser.add_argument("transcriptome")
    parser.add_argument("transcript_accession")
    args = parser.parse_args()

    bam_fn = Defaults.alignment_filename(args.workspace, args.sra_accession, args.transcriptome, "bam")
    sam_transcript_fn = Defaults.transcript_alignment_filename(args.workspace, args.sra_accession, args.transcriptome, args.transcript_accession, "sam")

    run_command("samtools", "view", bam_fn, args.transcript_accession, "-o", sam_transcript_fn)
    print(sam_transcript_fn)
