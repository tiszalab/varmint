"""
varmint CLI entry point.

This module exposes `main()` which is wired to the console script in pyproject.toml.
It calls `met_variant` and writes a TSV file.
"""
from __future__ import annotations

import argparse
import sys


MSG_PREFIX = "[varmint]"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="varmint",
        description=(
            "Compute allele frequencies from a BAM against a FASTA and annotate coding "
            "effects using a GFF (CDS). Outputs a TSV table."
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
        default=10,
        help="Minimum depth at a position to report variants",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()

    try:
        # Import at runtime so `varmint --help` works even if deps are missing
        from variant_funcs import met_variant  # type: ignore
    except Exception as e:
        sys.stderr.write(f"{MSG_PREFIX} ERROR importing variant_funcs: {e}\n")
        return 1

    try:
        df = met_variant(
            bam_path=args.bam,
            fasta_path=args.fasta,
            gff_path=args.gff,
            min_base_qual=args.min_base_qual,
            min_depth=args.min_depth,
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
