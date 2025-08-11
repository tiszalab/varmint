# function to assess variants from .bam
"""
Variant processing utilities for metagenomic alignments.

Provides `met_variant` which computes allele frequencies from a BAM
against a reference FASTA, and annotates coding consequences using a GFF
(CDS features) table. Returns a Polars DataFrame.
"""

from __future__ import annotations

from collections import Counter, defaultdict
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

import polars as pl
import pysam
from Bio import SeqIO
from Bio.Seq import Seq


# Complement lookup for single bases
_COMPLEMENT_TAB = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _comp_base(b: str) -> str:
    return b.translate(_COMPLEMENT_TAB)


def _lncomb(n: int, k: int) -> float:
    """Natural log of binomial coefficient C(n, k) using lgamma for stability."""
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def _fisher_exact_two_sided(a: int, b: int, c: int, d: int) -> float:
    """Two-sided Fisher's exact test p-value for 2x2 table [[a,b],[c,d]].

    Uses the standard definition: sum of probabilities of all tables with the same
    margins whose probability is less than or equal to the observed table's probability.
    Computed via hypergeometric probabilities in log-space for numerical stability.
    """
    r1 = a + b
    r2 = c + d
    c1 = a + c
    c2 = b + d
    n = r1 + r2
    # If any margin is zero, the table is degenerate; return p=1.0 (no evidence of bias)
    if r1 == 0 or r2 == 0 or c1 == 0 or c2 == 0:
        return 1.0
    lower = max(0, r1 - c2)
    upper = min(r1, c1)

    def logp(x: int) -> float:
        return _lncomb(c1, x) + _lncomb(c2, r1 - x) - _lncomb(n, r1)

    lp_obs = logp(a)
    # Precompute probabilities for all feasible tables
    lps = [logp(x) for x in range(lower, upper + 1)]
    # Compute two-sided p-value as sum of probs <= observed prob (with small tol)
    tol = 1e-12
    p = sum(math.exp(lp) for lp in lps if lp <= lp_obs + tol)
    # Guard against tiny numerical issues
    if p < 0:
        p = 0.0
    if p > 1:
        p = 1.0
    return p


@dataclass
class SegmentRec:
    contig: str
    transcript_id: str
    gene: Optional[str]
    start: int  # 1-based inclusive
    end: int  # 1-based inclusive
    strand: str  # '+' or '-'
    oriented_index: int  # index of this segment in transcriptional order
    cumulative_offset: int  # first base offset of this segment in CDS (0-based)
    length: int


@dataclass
class TranscriptModel:
    contig: str
    transcript_id: str
    gene: Optional[str]
    strand: str
    segments: List[SegmentRec]
    cds_seq: str  # 5'->3' coding DNA sequence (upper case)


def _parse_gff_attributes(attr_field: str) -> Dict[str, str]:
    """Parse GFF3 or GTF-like attributes field into a dict."""
    attrs: Dict[str, str] = {}
    s = attr_field.strip().strip(";")
    if not s:
        return attrs
    for item in s.split(";"):
        item = item.strip()
        if not item:
            continue
        if "=" in item:  # GFF3 style key=value
            k, v = item.split("=", 1)
            attrs[k.strip()] = v.strip().strip('"')
        else:  # GTF style: key "value"
            parts = item.split()
            if len(parts) >= 2:
                attrs[parts[0].strip()] = parts[1].strip().strip('"')
    return attrs


def _build_transcripts(
    gff_path: str, ref_seqs: Dict[str, str]
) -> Tuple[Dict[str, TranscriptModel], Dict[str, List[SegmentRec]]]:
    """Build transcript models (CDS only) and per-contig segment index.

    Returns (transcripts_by_id, segments_by_contig)
    """
    grouped: Dict[Tuple[str, str], List[Tuple[int, int, str, Optional[str]]]] = defaultdict(list)
    gene_name_by_key: Dict[Tuple[str, str], Optional[str]] = {}
    strand_by_key: Dict[Tuple[str, str], str] = {}

    with open(gff_path, "r") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            seqid, _source, ftype, start, end, _score, strand, _phase, attrs = parts
            if ftype != "CDS":
                continue
            start_i = int(start)
            end_i = int(end)
            attrs_d = _parse_gff_attributes(attrs)
            tid = (
                attrs_d.get("Parent")
                or attrs_d.get("transcript_id")
                or attrs_d.get("gene_id")
                or attrs_d.get("ID")
                or attrs_d.get("Name")
                or attrs_d.get("locus_tag")
            )
            if tid is None:
                tid = f"cds_{seqid}_{start}_{end}"
            gene_name = (
                attrs_d.get("gene")
                or attrs_d.get("Name")
                or attrs_d.get("locus_tag")
                or attrs_d.get("ID")
                or tid
            )
            key = (seqid, tid)
            grouped[key].append((start_i, end_i, strand, gene_name))
            strand_by_key[key] = strand
            gene_name_by_key[key] = gene_name

    transcripts: Dict[str, TranscriptModel] = {}
    segments_by_contig: Dict[str, List[SegmentRec]] = defaultdict(list)

    for (contig, tid), segs in grouped.items():
        if contig not in ref_seqs:
            # No sequence available for this contig in FASTA
            continue
        strand = strand_by_key[(contig, tid)]
        gene_name = gene_name_by_key[(contig, tid)]
        # Build coding sequence correctly:
        # - '+' strand: concatenate segments in ascending genomic order
        # - '-' strand: concatenate segments in ascending order then reverse-complement
        segs_asc = sorted(segs, key=lambda x: x[0])
        seq_pieces: List[str] = []
        for (s, e, _st, _gn) in segs_asc:
            s1 = max(1, s)
            e1 = e
            seq_pieces.append(ref_seqs[contig][s1 - 1 : e1])

        cds = "".join(seq_pieces)
        if strand == "-":
            cds = str(Seq(cds).reverse_complement())

        # Compute cumulative offsets in transcript 5'->3' order.
        # For '-' strand, transcript order corresponds to descending genomic order;
        # for '+', it corresponds to ascending genomic order.
        segs_for_offsets = (
            sorted(segs, key=lambda x: x[0]) if strand == "+" else sorted(segs, key=lambda x: x[0], reverse=True)
        )
        cum = 0
        seg_recs: List[SegmentRec] = []
        for idx, (s, e, _st, _gn) in enumerate(segs_for_offsets):
            seg_rec = SegmentRec(
                contig=contig,
                transcript_id=tid,
                gene=gene_name,
                start=s,
                end=e,
                strand=strand,
                oriented_index=idx,
                cumulative_offset=cum,
                length=e - s + 1,
            )
            cum += e - s + 1
            segments_by_contig[contig].append(seg_rec)
            seg_recs.append(seg_rec)

        tm = TranscriptModel(
            contig=contig,
            transcript_id=tid,
            gene=gene_name,
            strand=strand,
            segments=seg_recs,
            cds_seq=cds.upper(),
        )
        transcripts[tid] = tm

    # Sort segments on each contig by start for quicker linear scans
    for contig in list(segments_by_contig.keys()):
        segments_by_contig[contig].sort(key=lambda r: r.start)

    return transcripts, segments_by_contig


def _annotate_coding_effect(
    contig: str,
    pos1: int,
    alt_base: str,
    transcripts: Dict[str, TranscriptModel],
    segments_for_contig: List[SegmentRec],
) -> List[Dict[str, object]]:
    """Annotate a SNV at contig:pos1 (1-based) against overlapping CDS.

    Returns a list of annotations (one per overlapping transcript). Empty list if intergenic.
    """
    out: List[Dict[str, object]] = []
    for rec in segments_for_contig:
        if rec.start <= pos1 <= rec.end:
            tm = transcripts.get(rec.transcript_id)
            if tm is None or not tm.cds_seq:
                continue

            if tm.strand == "+":
                offset_in_seg = pos1 - rec.start
                alt_coding_base = alt_base
            else:
                offset_in_seg = rec.end - pos1
                alt_coding_base = _comp_base(alt_base)

            idx = rec.cumulative_offset + offset_in_seg
            if idx < 0 or idx >= len(tm.cds_seq):
                continue
            codon_start = (idx // 3) * 3
            if codon_start + 3 > len(tm.cds_seq):
                continue
            codon_ref = tm.cds_seq[codon_start : codon_start + 3]
            codon_list = list(codon_ref)
            codon_list[idx % 3] = alt_coding_base
            codon_alt = "".join(codon_list)

            # Translate using the Standard genetic code (NCBI table 1)
            aa_ref = str(Seq(codon_ref).translate(table=1))
            aa_alt = str(Seq(codon_alt).translate(table=1))

            if aa_ref == aa_alt:
                effect = "synonymous"
            elif aa_alt == "*":
                effect = "nonsense"
            elif aa_ref == "*":
                effect = "stop_loss"
            else:
                effect = "missense"

            out.append(
                {
                    "is_coding": True,
                    "gene": tm.gene,
                    "transcript_id": tm.transcript_id,
                    "strand": tm.strand,
                    "codon_ref": codon_ref,
                    "codon_alt": codon_alt,
                    "aa_ref": aa_ref,
                    "aa_alt": aa_alt,
                    "codon_index": (idx // 3) + 1,  # 1-based
                    "codon_pos": (idx % 3) + 1,  # 1..3
                    "effect": effect,
                }
            )

    return out


def met_variant(
    bam_path: str,
    fasta_path: str,
    gff_path: str,
    min_base_qual: int = 20,
    min_depth: int = 10,
) -> pl.DataFrame:
    """Compute allele frequencies per covered position and annotate coding effects.

    Parameters
    ----------
    bam_path : str
        Path to coordinate-sorted, indexed BAM file (.bam + .bai).
    fasta_path : str
        Reference genome FASTA used for alignment.
    gff_path : str
        Gene feature table (GFF3/GTF-like). CDS features are used for annotation.
    min_base_qual : int, optional
        Minimum base quality (Phred) to include a read base in counts. Default 20.
    min_depth : int, optional
        Minimum depth at a position to report allele frequencies. Default 10.

    Returns
    -------
    polars.DataFrame
        One row per detected alt allele (SNV) meeting thresholds, with allele
        frequency and coding consequence annotations when applicable.
        Includes per-allele average base qualities (avgq_A, avgq_C, avgq_G, avgq_T)
        at each reported position, computed over bases that pass filtering.
        A row is also emitted for positions with sufficient coverage even when
        no alternate allele is observed (reference-only sites; alt=None, af=0).
    """

    # Load reference sequences (upper-case strings)
    ref_seqs: Dict[str, str] = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")}

    # Build CDS transcript models and a per-contig segment index
    transcripts, segments_by_contig = _build_transcripts(gff_path, ref_seqs)

    # Open BAM
    bam = pysam.AlignmentFile(bam_path, "rb")

    records: List[Dict[str, object]] = []

    target_contigs = [c for c in bam.references if c in ref_seqs]
    for contig in target_contigs:
        # Iterate covered positions using pileup
        for puc in bam.pileup(
            contig,
            0,
            len(ref_seqs[contig]),
            truncate=True,
            stepper="samtools",
            min_base_quality=min_base_qual,
        ):
            pos1 = puc.reference_pos + 1  # convert to 1-based
            if pos1 < 1 or pos1 > len(ref_seqs[contig]):
                continue
            ref_base = ref_seqs[contig][pos1 - 1]
            if ref_base not in {"A", "C", "G", "T"}:
                # Skip ambiguous reference bases
                continue

            counts: Counter = Counter()
            base_fwd: Dict[str, int] = defaultdict(int)
            base_rev: Dict[str, int] = defaultdict(int)
            # Track summed qualities and counts per base to compute averages
            qual_sums: Dict[str, int] = defaultdict(int)
            qual_ns: Dict[str, int] = defaultdict(int)
            # Count bases meeting quality and not del/refskip
            for pr in puc.pileups:
                if pr.is_del or pr.is_refskip:
                    continue
                qp = pr.query_position
                if qp is None:
                    continue
                seq = pr.alignment.query_sequence
                if seq is None:
                    continue
                base = seq[qp]
                # Extra guard on base quality (min_base_quality already used in pileup)
                quals = pr.alignment.query_qualities
                if quals is not None and quals[qp] < min_base_qual:
                    continue
                bU = base.upper()
                if bU in {"A", "C", "G", "T"}:
                    counts[bU] += 1
                    if pr.alignment.is_reverse:
                        base_rev[bU] += 1
                    else:
                        base_fwd[bU] += 1
                    if quals is not None:
                        # Accumulate qualities for average calculations
                        qual_sums[bU] += int(quals[qp])
                        qual_ns[bU] += 1

            depth = sum(counts.get(b, 0) for b in ("A", "C", "G", "T"))
            if depth < min_depth:
                continue

            # Compute average qualities per nucleotide (None if no qualifying observations)
            avgq_A = (qual_sums["A"] / qual_ns["A"]) if qual_ns["A"] > 0 else None
            avgq_C = (qual_sums["C"] / qual_ns["C"]) if qual_ns["C"] > 0 else None
            avgq_G = (qual_sums["G"] / qual_ns["G"]) if qual_ns["G"] > 0 else None
            avgq_T = (qual_sums["T"] / qual_ns["T"]) if qual_ns["T"] > 0 else None

            # If no alt alleles observed (all reads agree with reference), emit one row
            has_alt = any(alt != ref_base and counts.get(alt, 0) > 0 for alt in ("A", "C", "G", "T"))
            if not has_alt:
                ref_fwd = base_fwd.get(ref_base, 0)
                ref_rev = base_rev.get(ref_base, 0)
                records.append(
                    {
                        "contig": contig,
                        "pos": pos1,
                        "ref": ref_base,
                        "alt": None,
                        "depth": depth,
                        "ref_count": counts.get(ref_base, 0),
                        "alt_count": 0,
                        "af": 0.0,
                        "count_A": counts.get("A", 0),
                        "count_C": counts.get("C", 0),
                        "count_G": counts.get("G", 0),
                        "count_T": counts.get("T", 0),
                        "avgq_A": avgq_A,
                        "avgq_C": avgq_C,
                        "avgq_G": avgq_G,
                        "avgq_T": avgq_T,
                        "ref_fwd": ref_fwd,
                        "ref_rev": ref_rev,
                        "alt_fwd": 0,
                        "alt_rev": 0,
                        "strand_bias_p": None,
                        "is_coding": False,
                        "gene": None,
                        "transcript_id": None,
                        "strand": None,
                        "codon_ref": None,
                        "codon_alt": None,
                        "aa_ref": None,
                        "aa_alt": None,
                        "codon_index": None,
                        "codon_pos": None,
                        "effect": None,
                    }
                )

            # Produce one row per alt allele observed
            for alt in ("A", "C", "G", "T"):
                c_alt = counts.get(alt, 0)
                if alt == ref_base or c_alt == 0:
                    continue
                af = c_alt / depth if depth > 0 else 0.0
                # Strand-specific counts and Fisher exact test
                alt_fwd = base_fwd.get(alt, 0)
                alt_rev = base_rev.get(alt, 0)
                ref_fwd = base_fwd.get(ref_base, 0)
                ref_rev = base_rev.get(ref_base, 0)
                sb_p = _fisher_exact_two_sided(alt_fwd, alt_rev, ref_fwd, ref_rev)

                # Coding consequence annotations (possibly multiple overlapping CDS)
                annot_list = _annotate_coding_effect(
                    contig,
                    pos1,
                    alt,
                    transcripts,
                    segments_by_contig.get(contig, []),
                )

                if not annot_list:
                    # Intergenic or non-CDS
                    records.append(
                        {
                            "contig": contig,
                            "pos": pos1,
                            "ref": ref_base,
                            "alt": alt,
                            "depth": depth,
                            "ref_count": counts.get(ref_base, 0),
                            "alt_count": c_alt,
                            "af": af,
                            "count_A": counts.get("A", 0),
                            "count_C": counts.get("C", 0),
                            "count_G": counts.get("G", 0),
                            "count_T": counts.get("T", 0),
                            "avgq_A": avgq_A,
                            "avgq_C": avgq_C,
                            "avgq_G": avgq_G,
                            "avgq_T": avgq_T,
                            "ref_fwd": ref_fwd,
                            "ref_rev": ref_rev,
                            "alt_fwd": alt_fwd,
                            "alt_rev": alt_rev,
                            "strand_bias_p": sb_p,
                            "is_coding": False,
                            "gene": None,
                            "transcript_id": None,
                            "strand": None,
                            "codon_ref": None,
                            "codon_alt": None,
                            "aa_ref": None,
                            "aa_alt": None,
                            "codon_index": None,
                            "codon_pos": None,
                            "effect": None,
                        }
                    )
                else:
                    for an in annot_list:
                        records.append(
                            {
                                "contig": contig,
                                "pos": pos1,
                                "ref": ref_base,
                                "alt": alt,
                                "depth": depth,
                                "ref_count": counts.get(ref_base, 0),
                                "alt_count": c_alt,
                                "af": af,
                                "count_A": counts.get("A", 0),
                                "count_C": counts.get("C", 0),
                                "count_G": counts.get("G", 0),
                                "count_T": counts.get("T", 0),
                                "avgq_A": avgq_A,
                                "avgq_C": avgq_C,
                                "avgq_G": avgq_G,
                                "avgq_T": avgq_T,
                                "ref_fwd": ref_fwd,
                                "ref_rev": ref_rev,
                                "alt_fwd": alt_fwd,
                                "alt_rev": alt_rev,
                                "strand_bias_p": sb_p,
                                "is_coding": an.get("is_coding", False),
                                "gene": an.get("gene"),
                                "transcript_id": an.get("transcript_id"),
                                "strand": an.get("strand"),
                                "codon_ref": an.get("codon_ref"),
                                "codon_alt": an.get("codon_alt"),
                                "aa_ref": an.get("aa_ref"),
                                "aa_alt": an.get("aa_alt"),
                                "codon_index": an.get("codon_index"),
                                "codon_pos": an.get("codon_pos"),
                                "effect": an.get("effect"),
                            }
                        )

    return pl.DataFrame(records)