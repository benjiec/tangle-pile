from pile import run_command, Defaults


def index(workspace, transcriptome, unclustered):
    if unclustered:
        fn = Defaults.transcriptome_fasta(workspace, transcriptome)
    else:
        fn = Defaults.transcriptome_cluster_fasta(workspace, transcriptome)

    run_command("bowtie2-build", fn, fn)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("workspace")
    parser.add_argument("transcriptome")
    parser.add_argument("--unclustered", action="store_true", default=False)
    args = parser.parse_args()

    index(args.workspace, args.transcriptome, args.unclustered)
