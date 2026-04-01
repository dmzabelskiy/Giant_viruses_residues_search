#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from _pipeline import project_root
from _family_pipeline import family_paths

PROJECT_ROOT = project_root()
INTERMEDIATE = PROJECT_ROOT / "data" / "intermediate"


def candidate_labels(rec: SeqRecord) -> list[str]:
    labels = []
    if rec.id:
        labels.append(str(rec.id).strip())
    desc = str(rec.description).strip()
    if desc and desc not in labels:
        labels.append(desc)
    return labels


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Extract viral and cellular members of a family from a supplied alignment.")
    ap.add_argument("--family", required=True, help="Exact gene_family value from curated_gene_sequences.tsv")
    ap.add_argument("--alignment", required=True, help="Path to the source family alignment FASTA")
    ap.add_argument("--results-subdir", default=None, help="Optional custom results directory name")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    family = args.family
    align_path = Path(args.alignment).expanduser().resolve()
    paths = family_paths(family, args.results_subdir)
    paths.ensure_dirs()

    if not align_path.exists():
        raise FileNotFoundError(f"Could not find alignment file: {align_path}")

    curated_all = pd.read_csv(INTERMEDIATE / "curated_gene_sequences.tsv", sep="\t")
    curated_viral = pd.read_csv(INTERMEDIATE / "curated_viral_sequences.tsv", sep="\t")
    curated_cell = pd.read_csv(INTERMEDIATE / "curated_cellular_sequences.tsv", sep="\t")

    fam_all = curated_all[curated_all["gene_family"] == family].copy()
    fam_viral = curated_viral[curated_viral["gene_family"] == family].copy()
    fam_cell = curated_cell[curated_cell["gene_family"] == family].copy()

    if fam_all.empty:
        raise ValueError(f"No curated rows found for family={family!r}")

    viral_labels = set(fam_viral["raw_label"].astype(str))
    cell_labels = set(fam_cell["raw_label"].astype(str))
    known_labels = set(fam_all["raw_label"].astype(str))

    records_all = list(SeqIO.parse(str(align_path), "fasta"))

    alignment_rows: list[dict[str, str]] = []
    viral_records_from_align: list[SeqRecord] = []
    cellular_records: list[SeqRecord] = []
    unmatched_records: list[SeqRecord] = []

    for rec in records_all:
        match_label = ""
        group = "unmatched"
        for label in candidate_labels(rec):
            if label in viral_labels:
                match_label = label
                group = "viral"
                break
            if label in cell_labels:
                match_label = label
                group = "cellular"
                break
            if label in known_labels:
                match_label = label
                group = "known_other"
                break

        rec2 = rec[:]
        if match_label:
            rec2.id = match_label
            rec2.name = match_label
        rec2.description = f"{rec2.id} source={group} family={family}"

        if group == "viral":
            viral_records_from_align.append(rec2)
        elif group == "cellular":
            cellular_records.append(rec2)
        else:
            unmatched_records.append(rec2)

        alignment_rows.append({
            "align_id": rec.id,
            "align_description": rec.description,
            "matched_label": match_label,
            "group": group,
        })

    out_map = paths.alignment_membership_tsv
    out_cell = paths.cellular_from_alignment_fasta
    out_viral_align = paths.viral_from_alignment_fasta
    out_unmatched = paths.alignment_unmatched_fasta
    out_summary = paths.alignment_summary_tsv

    pd.DataFrame(alignment_rows).to_csv(out_map, sep="\t", index=False)
    SeqIO.write(cellular_records, str(out_cell), "fasta")
    SeqIO.write(viral_records_from_align, str(out_viral_align), "fasta")
    SeqIO.write(unmatched_records, str(out_unmatched), "fasta")
    pd.DataFrame([
        {"metric": "family", "value": family},
        {"metric": "alignment_total", "value": len(records_all)},
        {"metric": "viral_from_alignment", "value": len(viral_records_from_align)},
        {"metric": "cellular_from_alignment", "value": len(cellular_records)},
        {"metric": "unmatched_or_other", "value": len(unmatched_records)},
    ]).to_csv(out_summary, sep="\t", index=False)

    print(f"Family: {family}")
    print(f"Results dir: {paths.results}")
    print(f"Intermediate dir: {paths.intermediate}")
    print(f"Alignment: {align_path}")
    print(f"Wrote: {out_map}")
    print(f"Wrote: {out_cell}")
    print(f"Wrote: {out_viral_align}")
    print(f"Wrote: {out_unmatched}")
    print(f"Wrote: {out_summary}")
    print("\nCounts:")
    print(f"  alignment total: {len(records_all)}")
    print(f"  viral from alignment: {len(viral_records_from_align)}")
    print(f"  cellular from alignment: {len(cellular_records)}")
    print(f"  unmatched/other from alignment: {len(unmatched_records)}")


if __name__ == "__main__":
    main()
