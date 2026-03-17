import os
import shutil
import tempfile
import subprocess
from pathlib import Path
from contextlib import contextmanager


def run_command(*cmd):
    subprocess.run(list(cmd), check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)

@contextmanager
def named_tempdir(dir=None, prefix=None, keep_on_error=False):
    with tempfile.TemporaryDirectory(dir=dir, prefix=prefix, delete=False) as temp_dir:
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



class Defaults(object):

    @contextmanager
    @staticmethod
    def named_tempdir(workspace):
        with named_tempdir(dir=workspace, prefix="temp-") as temp_dir:
            yield temp_dir

    @staticmethod
    def reads_dir(workspace):
        dn = f"{workspace}/reads"
        try:
            os.mkdir(dn)
        except FileExistsError as e:
            pass
        return dn

    @staticmethod
    def read_1(workspace, sra_accession):
        return f"{Defaults.reads_dir(workspace)}/{sra_accession}_1.fastq"

    @staticmethod
    def read_2(workspace, sra_accession):
        return f"{Defaults.reads_dir(workspace)}/{sra_accession}_2.fastq"
