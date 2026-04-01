#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from _pipeline import project_root
from _family_pipeline import family_paths

PROJECT_ROOT = project_root()
INTERMEDIATE = PROJECT_ROOT / "data" / "intermediate"
TOKEN_RE = re.compile(r"\b(contig_\d+_\d+)\b")


def parse_protein_token(raw_label: str) -> str:
    m = TOKEN_RE.search(str(raw_label).strip())
    return m.group(1) if m else ""


def load_fasta_records_by_token(fasta_path: Path, cache: dict[Path, dict[str, SeqRecord]]) -> dict[str, SeqRecord]:
    if fasta_path in cache:
        return cache[fasta_path]
    records: dict[str, SeqRecord] = {}
    for rec in SeqIO.parse(str(fasta_path), "fasta"):
        token = rec.description.split()[0]
        records[token] = rec
    cache[fasta_path] = records
    return records


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Extract viral sequences for an arbitrary gene family.")
    ap.add_argument("--family", required=True, help="Exact gene_family value from curated_gene_sequences.tsv")
    ap.add_argument("--results-subdir", default=None, help="Optional custom results directory name")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    family = args.family
    paths = family_paths(family, args.results_subdir)
    paths.ensure_dirs()

    curated = pd.read_csv(INTERMEDIATE / "curated_viral_sequences.tsv", sep="\t")
    mag_meta = pd.read_csv(INTERMEDIATE / "mag_metadata.tsv", sep="\t")

    sub = curated[curated["gene_family"] == family].copy()
    if sub.empty:
        known = sorted(curated["gene_family"].dropna().astype(str).unique())
        raise ValueError(f"No rows found for family={family!r}. Available examples: {known[:20]}")

    sub["protein_token"] = sub["raw_label"].map(parse_protein_token)
    fasta_lookup = (
        mag_meta[["mag_id_short", "fasta_path"]]
        .dropna()
        .drop_duplicates()
        .set_index("mag_id_short")["fasta_path"]
        .to_dict()
    )

    extracted_rows: list[dict[str, object]] = []
    extracted_records: list[SeqRecord] = []
    missing: list[dict[str, object]] = []
    cache: dict[Path, dict[str, SeqRecord]] = {}

    for _, row in sub.iterrows():
        mag_id = str(row.get("mag_id_short", "") or "")
        raw_label = str(row["raw_label"])
        protein_token = str(row["protein_token"])
        clade = str(row.get("clade", "") or "")
        fasta_path_str = fasta_lookup.get(mag_id, "")

        if not protein_token:
            missing.append({"reason": "failed_to_parse_protein_token", "mag_id_short": mag_id, "raw_label": raw_label})
            continue
        if not fasta_path_str:
            missing.append({"reason": "missing_fasta_for_mag", "mag_id_short": mag_id, "raw_label": raw_label, "protein_token": protein_token})
            continue

        fasta_path = Path(fasta_path_str)
        if not fasta_path.exists():
            missing.append({
                "reason": "fasta_file_not_found",
                "mag_id_short": mag_id,
                "raw_label": raw_label,
                "protein_token": protein_token,
                "fasta_path": str(fasta_path),
            })
            continue

        recs = load_fasta_records_by_token(fasta_path, cache)
        if protein_token not in recs:
            missing.append({
                "reason": "protein_token_not_found_in_fasta",
                "mag_id_short": mag_id,
                "raw_label": raw_label,
                "protein_token": protein_token,
                "fasta_file": fasta_path.name,
            })
            continue

        rec = recs[protein_token][:]
        clean_id = f"{mag_id}|{protein_token}|{clade}" if clade else f"{mag_id}|{protein_token}"
        rec.id = clean_id
        rec.name = clean_id
        rec.description = f"{clean_id} source=viral family={family} mag={mag_id} protein={protein_token} clade={clade}".strip()
        extracted_records.append(rec)
        extracted_rows.append({
            "clean_id": clean_id,
            "mag_id_short": mag_id,
            "protein_token": protein_token,
            "clade": clade,
            "raw_label": raw_label,
            "fasta_file": fasta_path.name,
            "length_aa": len(rec.seq),
        })

    out_faa = paths.viral_raw_fasta
    out_tsv = paths.viral_index_tsv
    out_missing = paths.viral_missing_tsv
    out_summary = paths.viral_extraction_summary_tsv

    SeqIO.write(extracted_records, str(out_faa), "fasta")
    pd.DataFrame(extracted_rows).to_csv(out_tsv, sep="\t", index=False)
    pd.DataFrame(missing).to_csv(out_missing, sep="\t", index=False)
    pd.DataFrame([
        {"metric": "family", "value": family},
        {"metric": "curated_viral_rows", "value": len(sub)},
        {"metric": "extracted_sequences", "value": len(extracted_records)},
        {"metric": "missing_rows", "value": len(missing)},
    ]).to_csv(out_summary, sep="\t", index=False)

    print(f"Family: {family}")
    print(f"Results dir: {paths.results}")
    print(f"Intermediate dir: {paths.intermediate}")
    print(f"Wrote: {out_faa}")
    print(f"Wrote: {out_tsv}")
    print(f"Wrote: {out_missing}")
    print(f"Wrote: {out_summary}")
    print(f"\nExtracted sequences: {len(extracted_records)}")
    print(f"Missing rows: {len(missing)}")
    if extracted_rows:
        df = pd.DataFrame(extracted_rows)
        print("\nClade counts:")
        print(df["clade"].replace("", "missing_clade").value_counts().to_string())


if __name__ == "__main__":
    main()
