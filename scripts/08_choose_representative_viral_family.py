#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
from Bio import SeqIO

from _family_pipeline import (
    RELAXED_CONSENSUS_CATEGORY,
    STRICT_CONSENSUS_CATEGORY,
    family_paths,
)

MISMATCH_PENALTY_FACTOR = 0.25
GAPLIKE_PENALTY_FACTOR = 0.50
GAPLIKE = set("-.X?")


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Choose a representative viral sequence for structure prediction.")
    ap.add_argument("--family", required=True)
    ap.add_argument("--results-subdir", default=None)
    return ap


def load_calls_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing residue-call table: {path}")
    df = pd.read_csv(path, sep="\t")
    required = {"alignment_pos_1based", "category", "top_viral_residue", "ranking_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Residue-call table is missing required columns: {sorted(missing)}")
    if "site_weight_score" not in df.columns:
        df["site_weight_score"] = df["ranking_score"]
    return df


def load_clean_viral_records(path: Path) -> Tuple[Dict[str, object], pd.DataFrame]:
    if not path.exists():
        raise FileNotFoundError(f"Missing viral FASTA: {path}")
    clean_records: Dict[str, object] = {}
    excluded_rows: List[dict] = []
    for rec in SeqIO.parse(str(path), "fasta"):
        seq = str(rec.seq)
        x_count = seq.count("X")
        if x_count > 0:
            excluded_rows.append({
                "seq_id": rec.id,
                "length": len(seq),
                "x_count": x_count,
                "x_fraction": x_count / len(seq) if seq else 0.0,
                "description": rec.description,
            })
            continue
        clean_records[rec.id] = rec
    excluded_df = pd.DataFrame(excluded_rows)
    return clean_records, excluded_df


def select_target_sites(calls: pd.DataFrame) -> pd.DataFrame:
    calls = calls.copy().sort_values("ranking_score", ascending=False)
    sel = calls[calls["category"].isin([STRICT_CONSENSUS_CATEGORY, RELAXED_CONSENSUS_CATEGORY])].copy()
    sel = sel[sel["top_viral_residue"].notna()].copy()
    sel["top_viral_residue"] = sel["top_viral_residue"].astype(str)
    sel = sel[sel["top_viral_residue"].str.len() == 1].copy()
    return sel


def build_site_targets(sel: pd.DataFrame) -> Tuple[Dict[int, str], Dict[int, float], Dict[int, str]]:
    target_aa: Dict[int, str] = {}
    weight: Dict[int, float] = {}
    category: Dict[int, str] = {}
    for _, row in sel.iterrows():
        pos0 = int(row["alignment_pos_1based"]) - 1
        aa = str(row["top_viral_residue"])
        w = float(row.get("site_weight_score", row["ranking_score"]))
        cat = str(row["category"])
        target_aa[pos0] = aa
        weight[pos0] = w
        category[pos0] = cat
    return target_aa, weight, category


def score_aligned_sequence(seq: str, target_aa: Dict[int, str], weight: Dict[int, float]) -> dict:
    matched = 0
    mismatched = 0
    gap_or_x = 0
    score = 0.0
    for pos0, expected_aa in target_aa.items():
        if pos0 >= len(seq):
            continue
        obs = seq[pos0]
        w = weight[pos0]
        if obs == expected_aa:
            matched += 1
            score += w
        elif obs in GAPLIKE:
            gap_or_x += 1
            score -= GAPLIKE_PENALTY_FACTOR * w
        else:
            mismatched += 1
            score -= MISMATCH_PENALTY_FACTOR * w
    return {"score": score, "matched_sites": matched, "mismatched_sites": mismatched, "gap_or_x_sites": gap_or_x}


def wrap_fasta_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))


def main() -> None:
    args = build_parser().parse_args()
    paths = family_paths(args.family, args.results_subdir)
    paths.ensure_dirs()
    align = paths.main_alignment_fasta if paths.main_alignment_fasta.exists() else paths.legacy_alignment_fasta
    viral_fasta = paths.viral_full_length_fasta
    calls_path = paths.residue_calls_tsv
    out_rank = paths.representative_ranking_tsv
    out_fasta = paths.representative_fasta
    out_excluded = paths.representative_excluded_tsv

    calls = load_calls_table(calls_path)
    clean_viral_raw, excluded_df = load_clean_viral_records(viral_fasta)

    if excluded_df.shape[0] > 0:
        excluded_df = excluded_df.sort_values(["x_count", "x_fraction", "length"], ascending=[False, False, False])
    excluded_df.to_csv(out_excluded, sep="\t", index=False)

    print(f"Viral full-length sequences without X: {len(clean_viral_raw)}")
    print(f"Excluded due to X: {len(excluded_df)}")
    print(f"Wrote: {out_excluded}")

    if not clean_viral_raw:
        raise ValueError("No clean viral full-length sequences remain after excluding X-containing records.")

    sel = select_target_sites(calls)
    if sel.empty:
        raise ValueError("No candidate viral-associated sites selected from the residue-call table.")

    target_aa, weight, site_category = build_site_targets(sel)

    aln_records = [rec for rec in SeqIO.parse(str(align), "fasta") if rec.id in set(clean_viral_raw.keys())]
    if not aln_records:
        raise ValueError("No clean viral sequences from viral full-length FASTA were found in the alignment.")

    rows: List[dict] = []
    for rec in aln_records:
        aligned_seq = str(rec.seq)
        raw_rec = clean_viral_raw[rec.id]
        raw_seq = str(raw_rec.seq).replace("-", "").replace(".", "")
        score_info = score_aligned_sequence(aligned_seq, target_aa, weight)

        high_match = moderate_match = high_total = moderate_total = 0
        for pos0, expected_aa in target_aa.items():
            if pos0 >= len(aligned_seq):
                continue
            obs = aligned_seq[pos0]
            cat = site_category[pos0]
            if cat == STRICT_CONSENSUS_CATEGORY:
                high_total += 1
                if obs == expected_aa:
                    high_match += 1
            elif cat == RELAXED_CONSENSUS_CATEGORY:
                moderate_total += 1
                if obs == expected_aa:
                    moderate_match += 1

        rows.append({
            "seq_id": rec.id,
            "score": score_info["score"],
            "matched_sites": score_info["matched_sites"],
            "mismatched_sites": score_info["mismatched_sites"],
            "gap_or_x_sites": score_info["gap_or_x_sites"],
            "high_confidence_matches": high_match,
            "high_confidence_total": high_total,
            "moderate_confidence_matches": moderate_match,
            "moderate_confidence_total": moderate_total,
            "high_confidence_match_fraction": high_match / high_total if high_total else 0.0,
            "moderate_confidence_match_fraction": moderate_match / moderate_total if moderate_total else 0.0,
            "raw_length": len(raw_seq),
            "description": raw_rec.description,
        })

    rank = pd.DataFrame(rows)
    rank = rank.sort_values(
        ["score", "matched_sites", "high_confidence_match_fraction", "moderate_confidence_match_fraction", "raw_length", "gap_or_x_sites"],
        ascending=[False, False, False, False, False, True],
    ).reset_index(drop=True)
    rank.to_csv(out_rank, sep="\t", index=False)

    best_id = rank.iloc[0]["seq_id"]
    best_rec = clean_viral_raw[best_id]
    best_seq = str(best_rec.seq).replace("-", "").replace(".", "")
    with open(out_fasta, "w") as fh:
        fh.write(f">{best_rec.id} viral_family_representative_for_colabfold\n")
        fh.write(wrap_fasta_sequence(best_seq) + "\n")

    print(f"\nWrote: {out_rank}")
    print(f"Wrote: {out_fasta}")
    print("\nTop 10 representative candidates:")
    print(rank[["seq_id", "score", "matched_sites", "mismatched_sites", "gap_or_x_sites", "high_confidence_matches", "high_confidence_total", "moderate_confidence_matches", "moderate_confidence_total", "raw_length"]].head(10).to_string(index=False))
    print(f"\nChosen representative: {best_id}")


if __name__ == "__main__":
    main()
