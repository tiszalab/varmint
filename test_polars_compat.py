from pathlib import Path

import polars as pl

from consensus_funcs import (
    _build_consensus_base_lookup,
    _get_majority_base,
    write_consensus_fasta,
)
from variant_funcs import _annotate_coding_effect, _build_transcripts


def _variant_frame() -> pl.DataFrame:
    return pl.DataFrame(
        [
            {
                "contig": "chr1",
                "pos": 1,
                "var_type": "REF",
                "allele_type": "ref",
                "ref_seq": "A",
                "alt_seq": None,
                "depth": 10,
                "AF": 0.4,
            },
            {
                "contig": "chr1",
                "pos": 1,
                "var_type": "SNV",
                "allele_type": "alt",
                "ref_seq": "A",
                "alt_seq": "G",
                "depth": 10,
                "AF": 0.6,
            },
            {
                "contig": "chr1",
                "pos": 2,
                "var_type": "REF",
                "allele_type": "ref",
                "ref_seq": "T",
                "alt_seq": None,
                "depth": 10,
                "AF": 1.0,
            },
        ]
    )


def test_consensus_dataframe_operations(tmp_path: Path) -> None:
    df = _variant_frame()

    assert _get_majority_base(df.filter(pl.col("pos") == 1), 0.5) == "G"
    assert _build_consensus_base_lookup(df, 0.5) == {
        ("chr1", 1): "G",
        ("chr1", 2): "T",
    }

    fasta_path = tmp_path / "ref.fasta"
    output_path = tmp_path / "consensus.fasta"
    fasta_path.write_text(">chr1\nATC\n")

    write_consensus_fasta(df, str(fasta_path), str(output_path))

    assert output_path.read_text().splitlines()[1] == "GTN"


def test_variant_annotation_and_dataframe_construction(tmp_path: Path) -> None:
    gff_path = tmp_path / "genes.gff"
    gff_path.write_text("chr1\ttest\tCDS\t1\t3\t.\t+\t0\tID=cds1\n")

    transcripts, segments = _build_transcripts(str(gff_path), {"chr1": "ATC"})
    annotations = _annotate_coding_effect(
        "chr1", 1, "G", transcripts, segments["chr1"]
    )
    frame = pl.DataFrame(annotations, infer_schema_length=len(annotations)).with_columns(
        pl.lit(None).cast(pl.Utf8).alias("VCF_PASS")
    )

    assert frame["effect"].to_list() == ["missense"]
    assert frame.schema["VCF_PASS"] == pl.String
