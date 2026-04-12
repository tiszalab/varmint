"""
varmint CLI entry point.

This module exposes `main()` which is wired to the console script in pyproject.toml.
writes TSV table.
"""
from __future__ import annotations

import argparse
import sys


MSG_PREFIX = "[varmint]"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="varmint",
        description=(
            "Variant/allele table generation and coding-effect annotation from BAM and VCF. "
            "Outputs TSV tables."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    p.add_argument("-b", "--bam", required=True, help="Input BAM file (coordinate-sorted, indexed)")
    p.add_argument(
        "-r",
        "--ref",
        "--fasta",
        dest="fasta",
        required=True,
        help="Reference FASTA used for alignment",
    )
    p.add_argument("-g", "--gff", required=True, help="GFF file with CDS features for annotation")
    p.add_argument("-o", "--out", required=True, help="Output TSV path")
    p.add_argument(
        "-v", 
        "--vcf", 
        type=str, 
        default=None,
        help="Input VCF/BCF with variants"
    )
    p.add_argument(
        "-q",
        "--min-base-qual",
        dest="min_base_qual",
        type=int,
        default=20,
        help="Minimum base quality (Phred) to count a base",
    )
    p.add_argument(
        "-d",
        "--min-depth",
        dest="min_depth",
        type=int,
        default=1,
        help="Minimum depth at a position to report variants",
    )
    p.add_argument(
        "--min-map-qual",
        dest="min_map_qual",
        type=int,
        default=0,
        help="Minimum mapping quality (MAPQ) to include a read",
    )
    p.add_argument(
        "--strand-bias-alpha",
        dest="strand_bias_alpha",
        type=float,
        default=None,
        help=(
            "If set, filter out alt alleles with Fisher's exact two-sided strand-bias p-value "
            "below this alpha (e.g., 0.001)."
        ),
    )
    p.add_argument(
        "-t",
        "--threads",
        dest="threads",
        type=int,
        default=1,
        help="Number of parallel threads to use for processing (default: 1)",
    )
    p.add_argument(
        "--consensus",
        dest="consensus",
        type=str,
        default=None,
        help="Optional path to write consensus FASTA file",
    )
    p.add_argument(
        "--consensus-af",
        dest="consensus_af",
        type=float,
        default=0.5,
        help="Allele frequency threshold for consensus IUPAC ambiguity codes (default: 0.5)",
    )

    return p.parse_args()

def main() -> int:
    args = parse_args()


    try:
        # Import at runtime so `varmint --help` works even if deps are missing
        from variant_funcs import met_variant_alleles
        from consensus_funcs import write_consensus_fasta
    except Exception as e:
        sys.stderr.write(f"{MSG_PREFIX} ERROR importing variant_funcs: {e}\n")
        return 1

    try:
        df = met_variant_alleles(
            bam_path=args.bam,
            fasta_path=args.fasta,
            gff_path=args.gff,
            vcf_path=args.vcf,
            min_base_qual=args.min_base_qual,
            min_depth=args.min_depth,
            min_map_qual=args.min_map_qual,
            strand_bias_alpha=args.strand_bias_alpha,
            threads=args.threads,
        )
    except Exception as e:
        sys.stderr.write(f"{MSG_PREFIX} ERROR during variant processing: {e}\n")
        return 1

    try:
        df.write_csv(args.out, separator="\t")
    except Exception as e:
        sys.stderr.write(f"{MSG_PREFIX} ERROR writing TSV to {args.out}: {e}\n")
        return 2

    sys.stdout.write(f"{MSG_PREFIX} Wrote {len(df)} rows to {args.out}\n")

    # Optional consensus FASTA output
    if args.consensus:
        try:
            write_consensus_fasta(
                df=df,
                fasta_path=args.fasta,
                output_path=args.consensus,
                af_threshold=args.consensus_af,
            )
            sys.stdout.write(f"{MSG_PREFIX} Wrote consensus FASTA to {args.consensus}\n")
        except Exception as e:
            sys.stderr.write(f"{MSG_PREFIX} ERROR writing consensus FASTA: {e}\n")
            return 3

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
