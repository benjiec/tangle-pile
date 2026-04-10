#!/usr/bin/env python3

from pathlib import Path
from pile import run_command
from pile.defaults import Defaults


def bowtie2_index(fn):
    run_command("bowtie2-build", fn, fn)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("fn")
    args = parser.parse_args()

    bowtie2_index(str(Path(Defaults.workspace_dir()) / args.fn))
