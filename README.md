# varmint

Metagenome variant processor: compute allele frequencies from BAM against FASTA and annotate coding effects using GFF CDS features.

## Install/run with uv

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
  --out variants.tsv \
  --min-base-qual 20 \
  --min-depth 10
```

## Dependencies

- polars (DataFrame/TSV)
- pysam (BAM pileups)
- biopython (FASTA/GFF parsing, translation)

Ensure your BAM is coordinate-sorted and indexed (`.bai`).
