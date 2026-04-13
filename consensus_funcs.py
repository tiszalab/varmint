# Consensus FASTA generation utilities for varmint
"""
Consensus FASTA generation from variant data.
"""

from typing import Dict, List, Optional, Tuple

import polars as pl
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


# IUPAC ambiguity codes for mixed bases
_IUPAC_AMBIGUITY = {
    frozenset(["A"]): "A",
    frozenset(["C"]): "C",
    frozenset(["G"]): "G",
    frozenset(["T"]): "T",
    frozenset(["A", "G"]): "R",
    frozenset(["C", "T"]): "Y",
    frozenset(["G", "C"]): "S",
    frozenset(["A", "T"]): "W",
    frozenset(["G", "T"]): "K",
    frozenset(["A", "C"]): "M",
    frozenset(["C", "G", "T"]): "B",
    frozenset(["A", "G", "T"]): "D",
    frozenset(["A", "C", "G"]): "V",
    frozenset(["A", "C", "G", "T"]): "N",
}


def _get_majority_base(pos_df: pl.DataFrame, af_threshold: float) -> str:
    """Determine consensus base for a position from allele data.
    
    Returns the majority base if AF >= threshold, otherwise returns
    IUPAC ambiguity code covering all alleles above (1 - threshold).
    """
    # Get reference row
    ref_row = pos_df.filter(pl.col("var_type") == "REF")
    if ref_row.is_empty():
        return "N"  # No coverage
    
    ref_base = ref_row["ref_seq"][0]
    depth = ref_row["depth"][0]
    
    if depth == 0:
        return "N"  # No coverage
    
    # Get all alt alleles and their frequencies
    alt_rows = pos_df.filter(pl.col("allele_type") == "alt")
    if alt_rows.is_empty():
        return ref_base  # No variants, return reference
    
    # Build frequency dict
    alleles = [(ref_base, ref_row["AF"][0])]
    for row in alt_rows.iter_rows(named=True):
        if row["alt_seq"] and len(row["alt_seq"]) == 1:  # SNV only
            alleles.append((row["alt_seq"], row["AF"]))
    
    # Find majority allele
    alleles.sort(key=lambda x: x[1], reverse=True)
    majority_base, majority_af = alleles[0]
    
    if majority_af >= af_threshold:
        return majority_base
    
    # Use ambiguity code: include all alleles > (1 - threshold)
    cutoff = 1.0 - af_threshold
    included = set()
    for base, af in alleles:
        if af > 0 and len(base) == 1:  # Only single-base alleles
            included.add(base.upper())
    
    if not included:
        return "N"
    
    return _IUPAC_AMBIGUITY.get(frozenset(included), "N")


def write_consensus_fasta(
    df: pl.DataFrame,
    fasta_path: str,
    output_path: str,
    af_threshold: float = 0.5,
) -> None:
    """Write consensus FASTA from variant DataFrame.
    
    Args:
        df: Polars DataFrame with variant data
        fasta_path: Path to reference FASTA (for sequence length and headers)
        output_path: Path to write consensus FASTA
        af_threshold: Minimum AF for unambiguous base call
    """
    # Load reference sequences
    ref_records = list(SeqIO.parse(fasta_path, "fasta"))
    ref_seqs: Dict[str, str] = {rec.id: str(rec.seq).upper() for rec in ref_records}
    
    # Build consensus records for all contigs
    consensus_records = []
    
    for ref_record in ref_records:
        contig_id = ref_record.id
        contig_df = df.filter(pl.col("contig") == contig_id)
        
        if contig_df.is_empty():
            # No data for this contig - write Ns
            consensus = "N" * len(ref_seqs[contig_id])
        else:
            # Build consensus position by position
            consensus_chars = []
            ref_len = len(ref_seqs[contig_id])
            
            # Get all positions that have data as a sorted list of Python ints
            positions_with_data = contig_df["pos"].unique().to_list()
            positions_set = set(int(p) for p in positions_with_data)
            
            # Create a list of (pos, df) pairs
            # Polars group_by yields (key_tuple, group_df) where key_tuple
            # is a tuple even for single-column groupby, e.g. (123,)
            pos_groups_list = []
            for key, group in contig_df.group_by("pos"):
                pos_int = int(key[0]) if isinstance(key, tuple) else int(key)
                pos_groups_list.append((pos_int, group))
            # Sort by position for efficient scanning
            pos_groups_list.sort(key=lambda x: x[0])
            
            # Build consensus sequence
            group_idx = 0
            for pos in range(1, ref_len + 1):
                # Fast path: check if position has data using Python int set
                if pos not in positions_set:
                    consensus_chars.append("N")
                    continue
                
                # Find the matching group (should be at or near current index)
                while group_idx < len(pos_groups_list) and pos_groups_list[group_idx][0] < pos:
                    group_idx += 1
                
                if group_idx < len(pos_groups_list) and pos_groups_list[group_idx][0] == pos:
                    _, pos_df = pos_groups_list[group_idx]
                    base = _get_majority_base(pos_df, af_threshold)
                    consensus_chars.append(base)
                else:
                    consensus_chars.append("N")  # No coverage
            
            consensus = "".join(consensus_chars)
        
        # Create record
        record = SeqRecord(
            Seq(consensus),
            id=contig_id,
            description=f"consensus af_threshold={af_threshold}",
        )
        consensus_records.append(record)
    
    # Write all contigs to output
    with open(output_path, "w") as fh:
        SeqIO.write(consensus_records, fh, "fasta")


def _build_consensus_base_lookup(
    df: pl.DataFrame,
    af_threshold: float,
) -> Dict[Tuple[str, int], str]:
    """Build a dict mapping (contig, pos) -> consensus base from variant DataFrame."""
    lookup: Dict[Tuple[str, int], str] = {}
    # Group by (contig, pos)
    for key, group in df.group_by(["contig", "pos"]):
        if isinstance(key, tuple):
            contig_val, pos_val = str(key[0]), int(key[1])
        else:
            continue
        base = _get_majority_base(group, af_threshold)
        lookup[(contig_val, pos_val)] = base
    return lookup


def write_consensus_aa_fasta(
    df: pl.DataFrame,
    fasta_path: str,
    gff_path: str,
    output_path: str,
    af_threshold: float = 0.5,
) -> None:
    """Write consensus amino acid FASTA for each coding region.

    For each transcript in the GFF, the consensus nucleotide sequence is
    assembled from the variant DataFrame at CDS positions, then translated.
    Positions without coverage are filled with N (producing X in translation).

    Args:
        df: Polars DataFrame with variant data (from met_variant_alleles)
        fasta_path: Path to reference FASTA
        gff_path: Path to GFF file with CDS features
        output_path: Path to write consensus amino acid FASTA
        af_threshold: Minimum AF for unambiguous base call
    """
    from variant_funcs import _build_transcripts

    # Load reference sequences
    ref_seqs: Dict[str, str] = {
        rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")
    }

    # Build transcript models from GFF
    transcripts, _segments_by_contig = _build_transcripts(gff_path, ref_seqs)

    if not transcripts:
        import sys
        sys.stderr.write(
            "[varmint] WARNING: No CDS features found in GFF. "
            "No consensus AA FASTA will be written.\n"
        )
        return

    # Build per-position consensus base lookup from the variant DataFrame
    base_lookup = _build_consensus_base_lookup(df, af_threshold)

    aa_records: List[SeqRecord] = []

    for composite_key, tm in sorted(transcripts.items(), key=lambda x: (x[0][0], x[0][1])):
        contig = tm.contig
        gene = tm.gene or composite_key[1]
        strand = tm.strand

        # Collect consensus bases at each CDS position in genomic order
        # Segments are stored in transcriptional order already
        cds_bases: List[str] = []
        for seg in tm.segments:
            if strand == "+":
                positions = range(seg.start, seg.end + 1)
            else:
                # For minus strand, segments are in descending genomic order
                # but we walk each segment ascending, then reverse-complement the whole CDS
                positions = range(seg.start, seg.end + 1)
            for gpos in positions:
                base = base_lookup.get((contig, gpos), "N")
                cds_bases.append(base)

        cds_nt = "".join(cds_bases)

        # For minus strand, reverse complement the assembled CDS
        if strand == "-":
            cds_nt = str(Seq(cds_nt).reverse_complement())

        # Translate (N codons -> X)
        # Pad to multiple of 3 if needed
        remainder = len(cds_nt) % 3
        if remainder != 0:
            cds_nt += "N" * (3 - remainder)

        aa_seq = str(Seq(cds_nt).translate())

        tid_str = composite_key[1]
        record = SeqRecord(
            Seq(aa_seq),
            id=f"{contig}|{tid_str}",
            description=f"gene={gene} strand={strand} af_threshold={af_threshold} nt_len={len(cds_nt)}",
        )
        aa_records.append(record)

    # Write all transcript AA sequences
    with open(output_path, "w") as fh:
        SeqIO.write(aa_records, fh, "fasta")
