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


def mkdir_exists(dn):
    if os.exists(dn):
        return
    try:
        os.mkdir(dn)
    except FileExistsError as e:
        pass


def clean_for_fn(fn):
    return re.sub("\W", "_", fn)


class Defaults(object):

    @staticmethod
    def workspace_dir(workspace):
        wsdir = os.environ.get('PILE_WORKSPACES_DIR')
        if wsdir is None or not wsdir.startswith("/"):
            raise Exception("PILE_WORKSPACES_DIR must be defined and must be an absolute path")
        if workspace.startswith("/"):
            raise Exception(f"Workspace {workspace} should not start with /")
        return str(Path(wsdir) / workspace)

    @contextmanager
    @staticmethod
    def named_tempdir(workspace):
        with named_tempdir(dir=Defaults.workspace_dir(workspace), prefix="temp-") as temp_dir:
            yield temp_dir

    @staticmethod
    def reads_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/reads"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def read_1(workspace, sra_accession):
        return f"{Defaults.reads_dir(workspace)}/{sra_accession}_1.fastq"

    @staticmethod
    def read_2(workspace, sra_accession):
        return f"{Defaults.reads_dir(workspace)}/{sra_accession}_2.fastq"

    @staticmethod
    def transcriptomes_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/transcriptomes"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def transcriptome_dir(workspace, transcriptome):
        clean_for_fn(transcriptome)
        dn = Path(Defaults.transcriptomes_dir(workspace)) / transcriptome
        mkdir_exists(dn)
        return dn

    @staticmethod
    def alignments_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/alignments"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def alignment_filename(workspace, sra_accession, transcriptome, suffix):
        dn = Defaults.alignments_dir(workspace)
        clean_for_fn(transcriptome)
        fn = f"{sra_accession}-{transcriptome}.{suffix}"
        return Path(dn) / fn

    @staticmethod
    def transcript_alignment_filename(workspace, sra_accession, transcriptome, transcript_accession, suffix):
        dn = Defaults.alignments_dir(workspace)
        clean_for_fn(transcriptome)
        clean_for_fn(transcript_accession)
        fn = f"{sra_accession}-{transcriptome}-{transcript_accession}.{suffix}"
        return Path(dn) / fn
