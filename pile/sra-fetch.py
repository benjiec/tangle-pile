from pile import run_command, named_tempdir, assert_exists, Defaults


def sra_fetch(workspace, sra_read_accession):
    with Defaults.named_tempdir(workspace) as tempdir:
        run_command("fasterq-dump", 
                    "-O", Defaults.reads_dir(workspace),
                    "-t", tempdir,
                    sra_read_accession)
        run_command("gzip", Defaults.read_1(workspace, sra_read_accession, gzip=False))
        run_command("gzip", Defaults.read_2(workspace, sra_read_accession, gzip=False))

    assert_exists(Defaults.read_1(workspace, sra_read_accession))
    assert_exists(Defaults.read_2(workspace, sra_read_accession))


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("sra_accession")
    args = parser.parse_args()

    sra_fetch(Defaults.workspace(), args.sra_read_accession)
