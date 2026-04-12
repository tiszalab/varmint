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
from concurrent.futures import ProcessPoolExecutor, as_completed

import polars as pl
import pysam
from Bio import SeqIO
from Bio.Seq import Seq


# Complement lookup for single bases
_COMPLEMENT_TAB = str.maketrans("ACGTNacgtn", "TGCANtgcan")

# Transition pairs (A<->G, C<->T)
_TRANSITION_PAIRS = frozenset([("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")])


def _comp_base(b: str) -> str:
    return b.translate(_COMPLEMENT_TAB)


def _classify_ts_tv(ref: str, alt: str) -> Optional[str]:
    """Classify SNV as transition ('Ts') or transversion ('Tv')."""
    if ref not in {"A", "C", "G", "T"} or alt not in {"A", "C", "G", "T"}:
        return None
    if (ref, alt) in _TRANSITION_PAIRS:
        return "Ts"
    return "Tv"


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
    transcript_id: Tuple[str, str]  # (contig, tid) composite key
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
    transcript_id: Tuple[str, str]  # (contig, tid) composite key
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
) -> Tuple[Dict[Tuple[str, str], TranscriptModel], Dict[str, List[SegmentRec]]]:
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

    transcripts: Dict[Tuple[str, str], TranscriptModel] = {}
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

        composite_key = (contig, tid)
        tm = TranscriptModel(
            contig=contig,
            transcript_id=composite_key,
            gene=gene_name,
            strand=strand,
            segments=seg_recs,
            cds_seq=cds.upper(),
        )
        transcripts[composite_key] = tm

    # Sort segments on each contig by start for quicker linear scans
    for contig in list(segments_by_contig.keys()):
        segments_by_contig[contig].sort(key=lambda r: r.start)

    return transcripts, segments_by_contig


def _annotate_coding_effect(
    contig: str,
    pos1: int,
    alt_base: str,
    transcripts: Dict[Tuple[str, str], TranscriptModel],
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


def _annotate_indel_effect(
    contig: str,
    pos1: int,
    ref_seq: str,
    alt_seq: str,
    transcripts: Dict[Tuple[str, str], TranscriptModel],
    segments_for_contig: List[SegmentRec],
) -> List[Dict[str, object]]:
    """Annotate an indel using VCF-style REF/ALT sequences anchored at pos1.

    Classification: frameshift vs in-frame insertion/deletion. Reports approximate
    codon index/position at the anchor base. Does not attempt full AA reconstruction.
    """
    out: List[Dict[str, object]] = []
    dlen = len(alt_seq) - len(ref_seq)
    for rec in segments_for_contig:
        if rec.start <= pos1 <= rec.end:
            tm = transcripts.get(rec.transcript_id)
            if tm is None or not tm.cds_seq:
                continue

            # Map anchor base into CDS index
            if tm.strand == "+":
                offset_in_seg = pos1 - rec.start
            else:
                offset_in_seg = rec.end - pos1
            idx = rec.cumulative_offset + offset_in_seg
            if idx < 0 or idx >= len(tm.cds_seq):
                continue

            if dlen % 3 != 0:
                effect = "frameshift"
            else:
                if dlen > 0:
                    effect = "inframe_insertion"
                elif dlen < 0:
                    effect = "inframe_deletion"
                else:
                    effect = "unknown"

            out.append(
                {
                    "is_coding": True,
                    "gene": tm.gene,
                    "transcript_id": tm.transcript_id,
                    "strand": tm.strand,
                    "codon_ref": None,
                    "codon_alt": None,
                    "aa_ref": None,
                    "aa_alt": None,
                    "codon_index": (idx // 3) + 1,
                    "codon_pos": (idx % 3) + 1,
                    "effect": effect,
                }
            )

    return out


def _process_region(
    contig: str,
    start: int,
    end: int,
    bam_path: str,
    ref_seq: str,
    ref_seqs: Dict[str, str],
    transcripts: Dict[Tuple[str, str], TranscriptModel],
    segments_for_contig: List[SegmentRec],
    vcf_filter_map: Dict[Tuple[str, int, str, str], str],
    min_base_qual: int,
    min_depth: int,
    min_map_qual: int,
    strand_bias_alpha: Optional[float],
) -> List[Dict[str, object]]:
    """Process a single genomic region and return variant records.

    This is a standalone top-level function to support ProcessPoolExecutor.
    """
    records: List[Dict[str, object]] = []

    def _vcf_lookup(ctg: str, pos1_l: int, ref_s: str, alt_s: str) -> Optional[str]:
        return vcf_filter_map.get((ctg, pos1_l, ref_s, alt_s))

    def _extract_insertion(read: pysam.AlignedSegment, ref_pos0: int, anchor_qpos: int, ins_len: int) -> Tuple[Optional[str], int, int]:
        """Return (inserted_seq, qual_sum, n_bases) for insertion after ref_pos0 in read.
        Uses CIGAR string to extract inserted bases directly from query_sequence.
        """
        if read.cigartuples is None or read.query_sequence is None:
            return None, 0, 0
        query_idx = 0
        ref_idx = read.reference_start
        for op, length in read.cigartuples:
            if ref_idx > ref_pos0:
                break
            if op == 0:  # M (match/mismatch)
                query_idx += length
                ref_idx += length
            elif op == 1:  # I (insertion)
                if ref_idx == ref_pos0 + 1 and anchor_qpos == query_idx - 1:
                    ins_start = query_idx
                    ins_end = min(query_idx + ins_len, len(read.query_sequence))
                    actual_len = ins_end - ins_start
                    if actual_len != ins_len:
                        return None, 0, 0
                    seq = read.query_sequence[ins_start:ins_end].upper()
                    qsum = 0
                    if read.query_qualities is not None:
                        for qpos in range(ins_start, ins_end):
                            qsum += int(read.query_qualities[qpos])
                    return seq, qsum, actual_len
                query_idx += length
            elif op == 2:  # D (deletion)
                ref_idx += length
            elif op == 4:  # S (soft clip)
                query_idx += length
            elif op == 3:  # N (skip)
                ref_idx += length
        return None, 0, 0

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for puc in bam.pileup(
            contig,
            start,
            end,
            truncate=True,
            stepper="samtools",
            min_base_quality=min_base_qual,
        ):
            pos1 = puc.reference_pos + 1
            if pos1 < 1 or pos1 > len(ref_seq):
                continue
            ref_base = ref_seq[pos1 - 1]
            if ref_base not in {"A", "C", "G", "T"}:
                continue

            # Accumulators
            base_counts: Counter = Counter()
            base_fwd: Dict[str, int] = defaultdict(int)
            base_rev: Dict[str, int] = defaultdict(int)
            bq_sum: Dict[str, int] = defaultdict(int)
            bq_n: Dict[str, int] = defaultdict(int)
            mq_sum: Dict[str, int] = defaultdict(int)
            mq_n: Dict[str, int] = defaultdict(int)

            ins_counts: Counter = Counter()
            ins_fwd: Dict[str, int] = defaultdict(int)
            ins_rev: Dict[str, int] = defaultdict(int)
            ins_qsum: Dict[str, int] = defaultdict(int)
            ins_qn: Dict[str, int] = defaultdict(int)
            ins_mqsum: Dict[str, int] = defaultdict(int)
            ins_mqn: Dict[str, int] = defaultdict(int)

            del_counts: Counter = Counter()
            del_fwd: Dict[Tuple[str, str], int] = defaultdict(int)
            del_rev: Dict[Tuple[str, str], int] = defaultdict(int)
            del_mqsum: Dict[Tuple[str, str], int] = defaultdict(int)
            del_mqn: Dict[Tuple[str, str], int] = defaultdict(int)

            depth = 0
            for pr in puc.pileups:
                if pr.is_refskip:
                    continue
                if pr.is_del:
                    continue
                read = pr.alignment
                if read.mapping_quality is not None and int(read.mapping_quality) < min_map_qual:
                    continue
                depth += 1

                # Deletion handling
                if pr.indel < 0:
                    k = -pr.indel
                    del_ref_seq = (ref_base + ref_seq[pos1 : pos1 + k]).upper()
                    del_alt_seq = ref_base
                    key = (del_ref_seq, del_alt_seq)
                    del_counts[key] += 1
                    if read.is_reverse:
                        del_rev[key] += 1
                    else:
                        del_fwd[key] += 1
                    del_mqsum[key] += int(read.mapping_quality)
                    del_mqn[key] += 1

                qp = pr.query_position
                if qp is None or read.query_sequence is None:
                    continue
                base = read.query_sequence[qp]
                quals = read.query_qualities
                if quals is not None and quals[qp] < min_base_qual:
                    continue
                bU = base.upper()
                if bU in {"A", "C", "G", "T"}:
                    base_counts[bU] += 1
                    if read.is_reverse:
                        base_rev[bU] += 1
                    else:
                        base_fwd[bU] += 1
                    if quals is not None:
                        bq_sum[bU] += int(quals[qp])
                        bq_n[bU] += 1
                    mq_sum[bU] += int(read.mapping_quality)
                    mq_n[bU] += 1

                # Insertion immediately after this base
                if pr.indel and pr.indel > 0:
                    ins_len = pr.indel
                    ins_seq, qsum, qn = _extract_insertion(read, puc.reference_pos, qp, ins_len)
                    if ins_seq and len(ins_seq) == ins_len:
                        ins_counts[ins_seq] += 1
                        if read.is_reverse:
                            ins_rev[ins_seq] += 1
                        else:
                            ins_fwd[ins_seq] += 1
                        if qn > 0:
                            ins_qsum[ins_seq] += qsum
                            ins_qn[ins_seq] += qn
                        ins_mqsum[ins_seq] += int(read.mapping_quality)
                        ins_mqn[ins_seq] += 1

            if depth < min_depth:
                continue

            ref_count = base_counts.get(ref_base, 0)
            ref_avgq = (bq_sum[ref_base] / bq_n[ref_base]) if bq_n[ref_base] > 0 else None
            ref_avgmq = (mq_sum[ref_base] / mq_n[ref_base]) if mq_n[ref_base] > 0 else None
            ref_fwd = base_fwd.get(ref_base, 0)
            ref_rev = base_rev.get(ref_base, 0)
            records.append(
                {
                    "contig": contig,
                    "pos": pos1,
                    "var_type": "REF",
                    "allele_type": "ref",
                    "ref_seq": ref_base,
                    "alt_seq": None,
                    "depth": depth,
                    "allele_count": ref_count,
                    "AF": ref_count / depth if depth > 0 else None,
                    "allele_avgq": ref_avgq,
                    "allele_avgmq": ref_avgmq,
                    "ref_fwd": ref_fwd,
                    "ref_rev": ref_rev,
                    "alt_fwd": None,
                    "alt_rev": None,
                    "strand_bias_p": None,
                    "ts_tv": None,
                    "VCF_PASS": _vcf_lookup(contig, pos1, ref_base, ref_base),
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

            # SNV alt rows
            for alt in ("A", "C", "G", "T"):
                if alt == ref_base:
                    continue
                c_alt = base_counts.get(alt, 0)
                if c_alt == 0:
                    continue
                alt_avgq = (bq_sum[alt] / bq_n[alt]) if bq_n[alt] > 0 else None
                alt_avgmq = (mq_sum[alt] / mq_n[alt]) if mq_n[alt] > 0 else None
                alt_fwd = base_fwd.get(alt, 0)
                alt_rev = base_rev.get(alt, 0)
                ref_fwd = base_fwd.get(ref_base, 0)
                ref_rev = base_rev.get(ref_base, 0)
                sb_p = _fisher_exact_two_sided(alt_fwd, alt_rev, ref_fwd, ref_rev)
                if strand_bias_alpha is not None and sb_p is not None and sb_p < strand_bias_alpha:
                    continue

                annot_list = _annotate_coding_effect(
                    contig, pos1, alt, transcripts, segments_for_contig,
                )
                vcf_status = _vcf_lookup(contig, pos1, ref_base, alt)
                if not annot_list:
                    records.append(
                        {
                            "contig": contig,
                            "pos": pos1,
                            "var_type": "SNV",
                            "allele_type": "alt",
                            "ref_seq": ref_base,
                            "alt_seq": alt,
                            "depth": depth,
                            "allele_count": c_alt,
                            "AF": c_alt / depth if depth > 0 else None,
                            "allele_avgq": alt_avgq,
                            "allele_avgmq": alt_avgmq,
                            "ref_fwd": ref_fwd,
                            "ref_rev": ref_rev,
                            "alt_fwd": alt_fwd,
                            "alt_rev": alt_rev,
                            "strand_bias_p": sb_p,
                            "ts_tv": _classify_ts_tv(ref_base, alt),
                            "VCF_PASS": vcf_status,
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
                                "var_type": "SNV",
                                "allele_type": "alt",
                                "ref_seq": ref_base,
                                "alt_seq": alt,
                                "depth": depth,
                                "allele_count": c_alt,
                                "AF": c_alt / depth if depth > 0 else None,
                                "allele_avgq": alt_avgq,
                                "allele_avgmq": alt_avgmq,
                                "ref_fwd": ref_fwd,
                                "ref_rev": ref_rev,
                                "alt_fwd": alt_fwd,
                                "alt_rev": alt_rev,
                                "strand_bias_p": sb_p,
                                "ts_tv": _classify_ts_tv(ref_base, alt),
                                "VCF_PASS": vcf_status,
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

            # Insertion alt rows
            for ins_seq, c_alt in ins_counts.items():
                alt_avgq = (ins_qsum[ins_seq] / ins_qn[ins_seq]) if ins_qn[ins_seq] > 0 else None
                alt_avgmq = (ins_mqsum[ins_seq] / ins_mqn[ins_seq]) if ins_mqn[ins_seq] > 0 else None
                alt_fwd = ins_fwd.get(ins_seq, 0)
                alt_rev = ins_rev.get(ins_seq, 0)
                ref_fwd = base_fwd.get(ref_base, 0)
                ref_rev = base_rev.get(ref_base, 0)
                sb_p = _fisher_exact_two_sided(alt_fwd, alt_rev, ref_fwd, ref_rev)
                if strand_bias_alpha is not None and sb_p is not None and sb_p < strand_bias_alpha:
                    continue

                ref_seq_for_ins = ref_base
                alt_full = ref_base + ins_seq
                vcf_status = _vcf_lookup(contig, pos1, ref_seq_for_ins, alt_full)
                annot_list = _annotate_indel_effect(
                    contig, pos1, ref_seq_for_ins, alt_full, transcripts, segments_for_contig,
                )
                if not annot_list:
                    records.append(
                        {
                            "contig": contig,
                            "pos": pos1,
                            "var_type": "INS",
                            "allele_type": "alt",
                            "ref_seq": ref_seq_for_ins,
                            "alt_seq": alt_full,
                            "depth": depth,
                            "allele_count": c_alt,
                            "AF": c_alt / depth if depth > 0 else None,
                            "allele_avgq": alt_avgq,
                            "allele_avgmq": alt_avgmq,
                            "ref_fwd": ref_fwd,
                            "ref_rev": ref_rev,
                            "alt_fwd": alt_fwd,
                            "alt_rev": alt_rev,
                            "strand_bias_p": sb_p,
                            "ts_tv": None,
                            "VCF_PASS": vcf_status,
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
                                "var_type": "INS",
                                "allele_type": "alt",
                                "ref_seq": ref_seq_for_ins,
                                "alt_seq": alt_full,
                                "depth": depth,
                                "allele_count": c_alt,
                                "AF": c_alt / depth if depth > 0 else None,
                                "allele_avgq": alt_avgq,
                                "allele_avgmq": alt_avgmq,
                                "ref_fwd": ref_fwd,
                                "ref_rev": ref_rev,
                                "alt_fwd": alt_fwd,
                                "alt_rev": alt_rev,
                                "strand_bias_p": sb_p,
                                "ts_tv": None,
                                "VCF_PASS": vcf_status,
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

            # Deletion alt rows
            for (del_ref_seq, del_alt_seq), c_alt in del_counts.items():
                alt_avgmq = (
                    del_mqsum.get((del_ref_seq, del_alt_seq), 0) / del_mqn.get((del_ref_seq, del_alt_seq), 0)
                    if del_mqn.get((del_ref_seq, del_alt_seq), 0) > 0
                    else None
                )
                alt_fwd = del_fwd.get((del_ref_seq, del_alt_seq), 0)
                alt_rev = del_rev.get((del_ref_seq, del_alt_seq), 0)
                ref_fwd = base_fwd.get(ref_base, 0)
                ref_rev = base_rev.get(ref_base, 0)
                sb_p = _fisher_exact_two_sided(alt_fwd, alt_rev, ref_fwd, ref_rev)
                if strand_bias_alpha is not None and sb_p is not None and sb_p < strand_bias_alpha:
                    continue

                vcf_status = _vcf_lookup(contig, pos1, del_ref_seq, del_alt_seq)
                annot_list = _annotate_indel_effect(
                    contig, pos1, del_ref_seq, del_alt_seq, transcripts, segments_for_contig,
                )
                if not annot_list:
                    records.append(
                        {
                            "contig": contig,
                            "pos": pos1,
                            "var_type": "DEL",
                            "allele_type": "alt",
                            "ref_seq": del_ref_seq,
                            "alt_seq": del_alt_seq,
                            "depth": depth,
                            "allele_count": c_alt,
                            "AF": c_alt / depth if depth > 0 else None,
                            "allele_avgq": None,
                            "allele_avgmq": alt_avgmq,
                            "ref_fwd": ref_fwd,
                            "ref_rev": ref_rev,
                            "alt_fwd": alt_fwd,
                            "alt_rev": alt_rev,
                            "strand_bias_p": sb_p,
                            "ts_tv": None,
                            "VCF_PASS": vcf_status,
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
                                "var_type": "DEL",
                                "allele_type": "alt",
                                "ref_seq": del_ref_seq,
                                "alt_seq": del_alt_seq,
                                "depth": depth,
                                "allele_count": c_alt,
                                "AF": c_alt / depth if depth > 0 else None,
                                "allele_avgq": None,
                                "allele_avgmq": alt_avgmq,
                                "ref_fwd": ref_fwd,
                                "ref_rev": ref_rev,
                                "alt_fwd": alt_fwd,
                                "alt_rev": alt_rev,
                                "strand_bias_p": sb_p,
                                "ts_tv": None,
                                "VCF_PASS": vcf_status,
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

    return records


def met_variant_alleles(
    bam_path: str,
    fasta_path: str,
    gff_path: str,
    vcf_path: Optional[str] = None,
    min_base_qual: int = 20,
    min_depth: int = 1,
    min_map_qual: int = 0,
    strand_bias_alpha: Optional[float] = None,
    threads: int = 1,
) -> pl.DataFrame:
    """One row per allele (including reference) per covered position.

    Columns: contig, pos, var_type (REF|SNV|INS|DEL), allele_type (ref|alt),
    ref_seq, alt_seq, depth, allele_count, allele_avgq, allele_avgmq, strand_bias_p, plus
    coding annotation fields (is_coding, gene, transcript_id, strand, codon_* , aa_*, effect).
    
    If strand_bias_alpha is provided, alt alleles (SNV/INS/DEL) with Fisher's exact
    two-sided p-value < strand_bias_alpha will be filtered out (not emitted).

    If vcf_path is provided, a column `VCF_PASS` will be added that reflects the
    FILTER field from the VCF for a matching allele at (contig, pos, ref, alt).
    'PASS' is recorded when no filters are present; otherwise the semicolon-joined
    list of filter names. Rows without a matching VCF allele will have VCF_PASS=None.
    """
    # Load reference sequences
    ref_seqs: Dict[str, str] = {rec.id: str(rec.seq).upper() for rec in SeqIO.parse(fasta_path, "fasta")}
    # Build CDS models
    transcripts, segments_by_contig = _build_transcripts(gff_path, ref_seqs)
    if not transcripts:
        import sys
        sys.stderr.write("[varmint] WARNING: No CDS features found in GFF. All positions will be non-coding.\n")

    # Optional: build VCF FILTER status map per allele
    vcf_filter_map: Dict[Tuple[str, int, str, str], str] = {}
    if vcf_path:
        try:
            # Suppress noisy htslib warnings (e.g., contig not defined in header)
            # while parsing VCF. This affects global htslib logging in this process.
            orig_verbosity = None
            try:
                orig_verbosity = pysam.get_verbosity()
                pysam.set_verbosity(0)
            except Exception:
                pass
            try:
                with pysam.VariantFile(vcf_path) as vcf:
                    for rec in vcf:  # iterate all records
                        cont = rec.contig
                        pos1_v = int(rec.pos)
                        ref_v = (rec.ref or "").upper()
                        alts_v = [a.upper() for a in (rec.alts or [])]
                        try:
                            filt_keys = list(rec.filter.keys())  # may be empty when PASS
                        except Exception:
                            filt_keys = []
                        status = "PASS" if not filt_keys or "PASS" in filt_keys else ";".join(sorted(str(x) for x in filt_keys))
                        for alt_v in alts_v:
                            vcf_filter_map[(cont, pos1_v, ref_v, alt_v)] = status
            finally:
                if orig_verbosity is not None:
                    try:
                        pysam.set_verbosity(orig_verbosity)
                    except Exception:
                        pass
        except Exception:
            # If VCF cannot be read, proceed without VCF_PASS annotations
            vcf_filter_map = {}

    # Determine target contigs
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        target_contigs = [c for c in bam.references if c in ref_seqs]

    # Build work chunks: (contig, start, end) tuples
    chunk_size = 10000  # Process in 10kb chunks for parallelism
    chunks = []
    for contig in target_contigs:
        contig_len = len(ref_seqs[contig])
        for start in range(0, contig_len, chunk_size):
            end = min(start + chunk_size, contig_len)
            chunks.append((contig, start, end))

    # Process chunks in parallel or serially
    records: List[Dict[str, object]] = []
    
    if threads > 1 and len(chunks) > 1:
        # Parallel processing with ProcessPoolExecutor
        # Prepare shared read-only data for pickling
        shared_data = {
            'bam_path': bam_path,
            'ref_seqs': ref_seqs,
            'transcripts': transcripts,
            'vcf_filter_map': vcf_filter_map,
            'min_base_qual': min_base_qual,
            'min_depth': min_depth,
            'min_map_qual': min_map_qual,
            'strand_bias_alpha': strand_bias_alpha,
        }
        
        with ProcessPoolExecutor(max_workers=threads) as executor:
            futures = {}
            for chunk in chunks:
                contig, start, end = chunk
                # Submit with contig-specific ref_seq and segments extracted
                future = executor.submit(
                    _process_region,
                    contig=contig,
                    start=start,
                    end=end,
                    ref_seq=ref_seqs[contig],
                    segments_for_contig=segments_by_contig.get(contig, []),
                    **shared_data
                )
                futures[future] = chunk
            for future in as_completed(futures):
                chunk_records = future.result()
                records.extend(chunk_records)
    else:
        # Serial processing
        for contig, start, end in chunks:
            chunk_records = _process_region(
                contig=contig,
                start=start,
                end=end,
                bam_path=bam_path,
                ref_seq=ref_seqs[contig],
                ref_seqs=ref_seqs,
                transcripts=transcripts,
                segments_for_contig=segments_by_contig.get(contig, []),
                vcf_filter_map=vcf_filter_map,
                min_base_qual=min_base_qual,
                min_depth=min_depth,
                min_map_qual=min_map_qual,
                strand_bias_alpha=strand_bias_alpha,
            )
            records.extend(chunk_records)
    
    # Sort records by contig and position for consistent output
    records.sort(key=lambda r: (r["contig"], r["pos"]))

    # Ensure consistent schema for optional VCF_PASS (may be Null for many rows).
    # Increase schema inference length so Polars sees string values and infers Utf8,
    # then explicitly cast VCF_PASS to Utf8 with nulls allowed.
    if records:
        df = pl.DataFrame(records, infer_schema_length=len(records))
    else:
        df = pl.DataFrame(records)
    if "VCF_PASS" in df.columns:
        df = df.with_columns(pl.col("VCF_PASS").cast(pl.Utf8))
    return df
