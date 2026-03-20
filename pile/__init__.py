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


class Defaults(object):

    @staticmethod
    def workspace():
        ws = os.environ.get('PILE_WORKSPACE')
        assert ws is not None
        return ws

    @staticmethod
    def workspace_dir(workspace=None):
        wsdir = os.environ.get('PILE_WORKSPACES_DIR')
        if wsdir is None or not wsdir.startswith("/"):
            raise Exception("PILE_WORKSPACES_DIR must be defined and must be an absolute path")
        if workspace is None:
            return wsdir
        if workspace.startswith("/"):
            raise Exception(f"Workspace {workspace} should not start with /")
        return str(Path(wsdir) / workspace)

    @contextmanager
    @staticmethod
    def named_tempdir(workspace):
        with named_tempdir(dir=Defaults.workspace_dir(workspace), prefix="temp-") as temp_dir:
            yield temp_dir

    @staticmethod
    def ncbi_dir():
        return Defaults.workspace_dir("ncbi")

    @staticmethod
    def ncbi_genome_genomic(genome_accession):
        return str(Path(Defaults.ncbi_dir()) / genome_accession / "genomic.fna")

    @staticmethod
    def ncbi_genome_proteins(genome_accession):
        return str(Path(Defaults.ncbi_dir()) / genome_accession / "proteins.faa")

    @staticmethod
    def reads_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/reads"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def read_1(workspace, sra_accession, gzip=True):
        fn = f"{Defaults.reads_dir(workspace)}/{sra_accession}_1.fastq"
        if gzip is False:
            return fn
        gzfn = Path(fn+".gz")
        if gzfn.exists():
            return str(gzfn)
        return fn

    @staticmethod
    def read_2(workspace, sra_accession, gzip=True):
        fn = f"{Defaults.reads_dir(workspace)}/{sra_accession}_2.fastq"
        if gzip is False:
            return fn
        gzfn = Path(fn+".gz")
        if gzfn.exists():
            return str(gzfn)
        return fn

    @staticmethod
    def transcriptomes_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/transcriptomes"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def transcriptome_dir(workspace, transcriptome):
        clean_for_fn(transcriptome)
        dn = str(Path(Defaults.transcriptomes_dir(workspace)) / transcriptome)
        mkdir_exists(dn)
        return dn

    @staticmethod
    def transcriptome_fasta(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "transcripts.fna")

    @staticmethod
    def transcriptome_proteins_fasta(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "proteins.faa")

    @staticmethod
    def transcriptome_cluster_fasta(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "transcript_clusters.fna")

    @staticmethod
    def transcriptome_filtered_read_1(workspace, transcriptome, sra_accession):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / f"{sra_accession}_filtered_1.fastq")

    @staticmethod
    def transcriptome_filtered_read_2(workspace, transcriptome, sra_accession):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / f"{sra_accession}_filtered_2.fastq")

    @staticmethod
    def alignments_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/alignments"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def alignment_filename(workspace, sra_accession, transcriptome, suffix):
        dn = Defaults.alignments_dir(workspace)
        clean_for_fn(transcriptome)
        fn = f"{sra_accession}.{transcriptome}.{suffix}"
        return str(Path(dn) / fn)

    @staticmethod
    def transcripts_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/transcripts"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def transcript_alignment_filename(workspace, sra_accession, transcriptome, transcript_accession, suffix):
        dn = Defaults.transcripts_dir(workspace)
        clean_for_fn(transcriptome)
        clean_for_fn(transcript_accession)
        fn = f"{transcript_accession}.{sra_accession}.{transcriptome}.{suffix}"
        return str(Path(dn) / fn)
