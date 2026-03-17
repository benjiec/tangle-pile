from pile import run_command, named_tempdir, assert_exists, Defaults


def sra_fetch(workspace, sra_read_accession):
    with Defaults.named_tempdir(workspace) as tempdir:
        run_command("/Users/benjie/sratoolkit.3.3.0-mac-arm64/bin/fasterq-dump.3", 
                    "-O", Defaults.reads_dir(workspace),
                    "-t", tempdir,
                    sra_read_accession)

    assert_exists(Defaults.read_1(workspace, sra_read_accession))
    assert_exists(Defaults.read_2(workspace, sra_read_accession))



if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("workspace")
    parser.add_argument("sra_read_accession")
    args = parser.parse_args()

    sra_fetch(args.workspace, args.sra_read_accession)
