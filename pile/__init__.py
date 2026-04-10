import re
import os
import sys
import shutil
import tempfile
import subprocess
from pathlib import Path
from contextlib import contextmanager


def run_command(*cmd):
    print(" ".join(list(cmd)))
    subprocess.run(list(cmd), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

@contextmanager
def named_tempdir(dir=None, prefix=None, keep_on_error=False):
    with tempfile.TemporaryDirectory(dir=dir, prefix=prefix) as temp_dir:
        try:
            yield temp_dir
        finally: 
            if not keep_on_error:
                shutil.rmtree(temp_dir)

@contextmanager
def named_tempfile(dir=None, suffix=None, keep_on_error=False):
    with tempfile.NamedTemporaryFile(dir=dir, suffix=suffix, delete=False) as tmpf:
        try:
            yield tmpf.name
        finally: 
            if not keep_on_error:
                os.remove(tmpf.name)

def assert_exists(fn):
    path = Path(fn)
    assert path.exists()

def mkdir_exists(dn):
    if os.path.exists(dn):
        return
    try:
        os.mkdir(dn)
    except FileExistsError as e:
        pass

def clean_for_fn(fn):
    return re.sub("\W", "_", fn)


def process_file_or_literal(conditional_for_literal, argument, task):
    if conditional_for_literal:
        task(argument)
    else:
        if argument == "-":
            input_file = sys.stdin
        else:
            input_file = open(argument, "r")
        for line in input_file:
            task(line.strip())


class FastaHeaderProvenance(object):

    @staticmethod
    def add_and_str(target_accession, provenance):
        assert "<" not in target_accession
        new_provenance = [target_accession]
        if provenance is not None:
            new_provenance.extend(provenance)
        return " < ".join(new_provenance)

    @staticmethod
    def unpack(header):
        if "<" not in header:
            return [header]
        tokens = re.split(r"\s*\<\s*", header)
        return tokens
