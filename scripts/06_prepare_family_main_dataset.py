#!/usr/bin/env python3
from __future__ import annotations

import argparse

from Bio import SeqIO

from _family_pipeline import family_paths


def ungap_record(rec):
    rec2 = rec[:]
    seq = str(rec.seq).replace("-", "").replace(".", "")
    rec2.seq = rec.seq.__class__(seq)
    return rec2


def assert_unique_ids(records, label: str) -> None:
    ids = [r.id for r in records]
    dup = sorted({x for x in ids if ids.count(x) > 1})
    if dup:
        raise ValueError(f"Duplicate IDs in {label}: {dup[:10]}")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Prepare the combined viral+cellular ungapped dataset for alignment.")
    ap.add_argument("--family", required=True)
    ap.add_argument("--results-subdir", default=None)
    return ap


def main() -> None:
    args = build_parser().parse_args()
    paths = family_paths(args.family, args.results_subdir)
    paths.ensure_dirs()

    viral_full_fasta = paths.viral_full_length_fasta
    cell_align_fasta = paths.cellular_from_alignment_fasta
    cell_ungapped_fasta = paths.cellular_ungapped_fasta
    combined_fasta = paths.main_ungapped_fasta

    viral_records = list(SeqIO.parse(str(viral_full_fasta), "fasta"))
    cell_records = [ungap_record(r) for r in SeqIO.parse(str(cell_align_fasta), "fasta")]
    assert_unique_ids(viral_records, "viral full-length FASTA")
    assert_unique_ids(cell_records, "cellular alignment FASTA")
    overlap = sorted(set(r.id for r in viral_records) & set(r.id for r in cell_records))
    if overlap:
        raise ValueError(f"Viral/cellular ID overlap detected: {overlap[:10]}")

    SeqIO.write(cell_records, str(cell_ungapped_fasta), "fasta")
    SeqIO.write(viral_records + cell_records, str(combined_fasta), "fasta")

    print(f"Family: {args.family}")
    print(f"Results dir: {paths.results}")
    print(f"Intermediate dir: {paths.intermediate}")
    print(f"Wrote: {cell_ungapped_fasta}")
    print(f"Wrote: {combined_fasta}")
    print("\nCounts:")
    print(f"  viral full-length-like: {len(viral_records)}")
    print(f"  cellular ungapped: {len(cell_records)}")
    print(f"  combined: {len(viral_records) + len(cell_records)}")


if __name__ == "__main__":
    main()
