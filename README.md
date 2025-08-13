# varmint

Sequence variant processor: compute allele frequencies from BAM against FASTA and annotate coding effects using GFF CDS features. Add a VCF to report on "statistically significant" variants from your favorite variant caller.

## Quick install with `pip`

*This will not install dependencies*

*This will install `varmint` into whichever environment (e.g. conda env) that you are in.*

Dependencies:
- polars (DataFrame/TSV)
- pysam (BAM pileups)
- biopython (FASTA/GFF parsing, translation)

```bash
cd /path/to/varmint/

pip install .
```

test that `varmint` is accessible now.

```bash
varmint -h
```

## Install/run with `uv`

Option A: Install as a tool (adds `varmint` to your PATH)

```bash
uv tool install .
varmint --help
```

Option B: Run within a project environment (no global tool install)

```bash
uv run python -m varmint_cli --help
uv run varmint -b sample.bam -r ref.fasta -g genes.gff -o variants.tsv
```

Option C: Editable dev install in a uv-managed venv

```bash
uv venv
source .venv/bin/activate
uv pip install -e .
varmint --help
```

## Usage

```bash
varmint \
  --bam sample.bam \
  --ref reference.fasta \
  --gff genes.gff \
  --vcf variants.vcf
  --out variant_stats.tsv \
  --min-base-qual 20 \
  --min-depth 10
```

**Ensure your BAM is coordinate-sorted and indexed (`.bai`).**

## Output Data Dictionary

| Column | Type | Description |
|---|---|---|
| contig | str | Reference contig/chromosome ID. Must match the FASTA. |
| pos | int (1-based) | Genomic position (1-based). |
| var_type | str | Variant type: REF (reference row), SNV, INS, DEL. |
| allele_type | str | Allele row type: ref (reference allele) or alt (variant allele). |
| ref_seq | str | Reference base for SNVs, or left-anchored reference allele for indels (e.g., A for SNV; ACG for deletion of CG). |
| alt_seq | str or null | Alt base for SNVs, or left-anchored alt allele for indels (e.g., AT for insertion of T after A). Null for reference rows. |
| depth | int or null | Total read depth at the position after filtering (MAPQ, base qual, etc.). |
| allele_count | int or null | Number of reads supporting this allele (for REF rows, count of the reference base). |
| allele_avgq | float or null | Mean base quality of bases supporting this allele (available for SNVs and insertions; deletions are null). |
| allele_avgmq | float or null | Mean mapping quality of reads supporting this allele. |
| strand_bias_p | float or null | Fisher’s exact two-sided p-value for strand bias of the alt vs ref (null for ref rows or insufficient data). |
| VCF_PASS | str or null | From input VCF FILTER for matching (contig, pos, ref, alt): "PASS" or semicolon-joined filter names. Null if no matching VCF allele or when no VCF provided. |
| is_coding | bool | True if allele overlaps a CDS feature. |
| gene | str or null | Gene name/ID for overlapping CDS (if available). |
| transcript_id | str or null | Transcript ID for overlapping CDS. |
| strand | "+"/"-" or null | Strand of the overlapping CDS. |
| codon_ref | str or null | Reference codon (for coding SNVs/indels when determinable). |
| codon_alt | str or null | Alternate codon (for coding SNVs/indels when determinable). |
| aa_ref | str or null | Reference amino acid (single-letter; when determinable). |
| aa_alt | str or null | Alternate amino acid (single-letter; when determinable). |
| codon_index | int or null | 1-based index of the affected codon within the CDS segment. |
| codon_pos | int (1–3) or null | Position within the codon (1, 2, or 3). |
| effect | str or null | Predicted effect (e.g., synonymous, missense, nonsense, frameshift, inframe_ins/del, mnp, unknown). |

