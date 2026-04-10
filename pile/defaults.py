import os
from pathlib import Path
from contextlib import contextmanager
from tangle.defaults import PathDefaultsBase


class Defaults(PathDefaultsBase):

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
    def transcriptome_fasta(workspace, transcriptome, gzip=True):
        fn = str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "transcripts.fna")
        if gzip is False:
            return fn
        gzfn = Path(fn+".gz")
        if gzfn.exists():
            return str(gzfn)
        return fn

    @staticmethod
    def transcriptome_salmon_index(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "transcripts.fna.salmon_index")

    @staticmethod
    def transcriptome_proteins_fasta(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "proteins.faa")

    @staticmethod
    def transcriptome_proteins_gff(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "proteins.gff3")

    @staticmethod
    def transcriptome_unclustered_fasta(workspace, transcriptome, gzip=True):
        fn = str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "transcripts.unclustered.fna")
        if gzip is False:
            return fn
        gzfn = Path(fn+".gz")
        if gzfn.exists():
            return str(gzfn)
        return fn

    @staticmethod
    def transcriptome_unclustered_salmon_index(workspace, transcriptome):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / "transcripts.unclustered.fna.salmon_index")

    @staticmethod
    def transcriptome_filtered_read_1(workspace, transcriptome, sra_accession):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / f"{sra_accession}_filtered_1.fastq")

    @staticmethod
    def transcriptome_filtered_read_2(workspace, transcriptome, sra_accession):
        return str(Path(Defaults.transcriptome_dir(workspace, transcriptome)) / f"{sra_accession}_filtered_2.fastq")

    @staticmethod
    def quants_dir(workspace):
        dn = f"{Defaults.workspace_dir(workspace)}/quants"
        mkdir_exists(dn)
        return dn

    @staticmethod
    def quant_transcriptome_dir(workspace, transcriptome):
        dn = str(Path(Defaults.quants_dir(workspace)) / transcriptome)
        mkdir_exists(dn)
        return dn

    @staticmethod
    def quant_dir(workspace, transcriptome, sra_accession):
        dn = str(Path(Defaults.quant_transcriptome_dir(workspace, transcriptome)) / sra_accession)
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


if __name__ == "__main__":
    Defaults.main(Defaults)
