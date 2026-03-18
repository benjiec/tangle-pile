from pile import run_command, Defaults, named_tempfile


def bowtie2_align(workspace, transcriptome, unclustered, sra_accession, sam_fn, cpus):
    reads_1 = Defaults.read_1(workspace, sra_accession)
    reads_2 = Defaults.read_2(workspace, sra_accession)
    if unclustered:
        tx = Defaults.transcriptome_fasta(workspace, transcriptome)
    else:
        tx = Defaults.transcriptome_cluster_fasta(workspace, transcriptome)

    run_command(
      "bowtie2", "--local", "--very-sensitive-local", "-p", str(cpus),
      "-x", tx,
      "-1", reads_1,
      "-2", reads_2,
      "-S", sam_fn
    )


def sam_to_sorted_bam(temp_dir, sam_fn, bam_fn):
    with named_tempfile(dir=temp_dir) as unsorted_bam_fn:
        run_command("samtools", "view", "-S", "-b", "-o", unsorted_bam_fn, sam_fn)
        run_command("samtools", "sort", unsorted_bam_fn, "-o", bam_fn)
        run_command("samtools", "index", bam_fn)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("sra_accession")
    parser.add_argument("transcriptome")
    parser.add_argument("--unclustered", action="store_true", default=False)
    parser.add_argument("--cpus", default=2)
    args = parser.parse_args()

    bam_fn = Defaults.alignment_filename(Defaults.workspace(), args.sra_accession, args.transcriptome, "bam")

    with Defaults.named_tempdir(Defaults.workspace()) as temp_dir:
        with named_tempfile(dir=temp_dir) as sam_fn:
            bowtie2_align(Defaults.workspace(), args.transcriptome, args.unclustered, args.sra_accession, sam_fn, args.cpus)
            sam_to_sorted_bam(temp_dir, sam_fn, bam_fn)
