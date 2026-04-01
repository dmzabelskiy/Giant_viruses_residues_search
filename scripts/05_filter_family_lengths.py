#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from _family_pipeline import family_paths



def ungapped_length(seq: str) -> int:
    return len(seq.replace("-", "").replace(".", ""))


def compute_auto_min_length(cell_lengths: list[int], quantile: float = 0.10, mad_scale: float = 3.0) -> tuple[int, dict[str, float]]:
    if not cell_lengths:
        raise ValueError("No cellular lengths available for automatic thresholding.")
    s = pd.Series(cell_lengths, dtype=float)
    median = float(s.median())
    mad = float((s - median).abs().median())
    q = float(s.quantile(quantile))
    robust_floor = median - mad_scale * mad
    chosen = int(math.floor(max(q, robust_floor, 1.0)))
    stats = {
        "n": int(s.shape[0]),
        "min": float(s.min()),
        "q05": float(s.quantile(0.05)),
        "q10": float(s.quantile(0.10)),
        "q25": float(s.quantile(0.25)),
        "median": median,
        "q75": float(s.quantile(0.75)),
        "q90": float(s.quantile(0.90)),
        "max": float(s.max()),
        "mad": mad,
        "robust_floor_median_minus_3mad": robust_floor,
        "chosen_min_length": float(chosen),
    }
    return chosen, stats


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Filter viral sequences by family-specific full-length threshold.")
    ap.add_argument("--family", required=True, help="Exact gene_family value")
    ap.add_argument("--results-subdir", default=None, help="Optional custom results directory name")
    ap.add_argument("--min-length", type=int, default=None, help="Manual override. If omitted, infer from cellular lengths.")
    ap.add_argument("--quantile", type=float, default=0.10, help="Lower cellular quantile used by auto-thresholding.")
    ap.add_argument("--mad-scale", type=float, default=3.0, help="Use median - mad_scale*MAD as a robust lower floor.")
    return ap


def main() -> None:
    args = build_parser().parse_args()
    paths = family_paths(args.family, args.results_subdir)
    paths.ensure_dirs()

    viral_index = paths.viral_index_tsv
    viral_fasta = paths.viral_raw_fasta
    cell_align_fasta = paths.cellular_from_alignment_fasta

    if not viral_index.exists():
        raise FileNotFoundError(f"Missing viral index: {viral_index}")
    if not viral_fasta.exists():
        raise FileNotFoundError(f"Missing viral FASTA: {viral_fasta}")
    if not cell_align_fasta.exists():
        raise FileNotFoundError(f"Missing cellular alignment FASTA: {cell_align_fasta}")

    idx = pd.read_csv(viral_index, sep="\t")
    cell_lengths = [ungapped_length(str(rec.seq)) for rec in SeqIO.parse(str(cell_align_fasta), "fasta")]
    if args.min_length is None:
        min_len, stats = compute_auto_min_length(cell_lengths, quantile=args.quantile, mad_scale=args.mad_scale)
        method = f"auto:max(q{int(args.quantile*100):02d}, median-{args.mad_scale}*MAD)"
    else:
        min_len = int(args.min_length)
        _, stats = compute_auto_min_length(cell_lengths, quantile=args.quantile, mad_scale=args.mad_scale)
        stats["chosen_min_length"] = float(min_len)
        method = "manual"

    full_ids = set(idx.loc[idx["length_aa"] >= min_len, "clean_id"])
    frag_ids = set(idx.loc[idx["length_aa"] < min_len, "clean_id"])
    full_records = []
    frag_records = []
    for rec in SeqIO.parse(str(viral_fasta), "fasta"):
        if rec.id in full_ids:
            full_records.append(rec)
        elif rec.id in frag_ids:
            frag_records.append(rec)

    full_df = idx[idx["length_aa"] >= min_len].copy()
    frag_df = idx[idx["length_aa"] < min_len].copy()

    full_fasta = paths.viral_full_length_fasta
    frag_fasta = paths.viral_fragments_fasta
    full_index = paths.viral_full_length_index_tsv
    frag_index = paths.viral_fragments_tsv
    summary_out = paths.viral_length_filter_summary_tsv

    SeqIO.write(full_records, str(full_fasta), "fasta")
    SeqIO.write(frag_records, str(frag_fasta), "fasta")
    full_df.to_csv(full_index, sep="\t", index=False)
    frag_df.to_csv(frag_index, sep="\t", index=False)

    summary_rows = [
        {"metric": "family", "value": args.family},
        {"metric": "threshold_method", "value": method},
        {"metric": "min_length_threshold", "value": min_len},
        {"metric": "viral_total", "value": len(idx)},
        {"metric": "viral_full_length_like", "value": len(full_df)},
        {"metric": "viral_fragments", "value": len(frag_df)},
    ]
    summary_rows.extend({"metric": f"cellular_length_{k}", "value": v} for k, v in stats.items())
    pd.DataFrame(summary_rows).to_csv(summary_out, sep="\t", index=False)

    print(f"Family: {args.family}")
    print(f"Results dir: {paths.results}")
    print(f"Intermediate dir: {paths.intermediate}")
    print(f"Threshold method: {method}")
    print(f"Minimum length threshold: {min_len} aa")
    print(f"Wrote: {full_fasta}")
    print(f"Wrote: {frag_fasta}")
    print(f"Wrote: {full_index}")
    print(f"Wrote: {frag_index}")
    print(f"Wrote: {summary_out}")
    print("\nCellular length stats:")
    for k in ["n", "min", "q05", "q10", "q25", "median", "q75", "q90", "max", "mad", "robust_floor_median_minus_3mad"]:
        print(f"  {k}: {stats[k]}")
    print("\nCounts:")
    print(f"  total viral: {len(idx)}")
    print(f"  full-length-like: {len(full_df)}")
    print(f"  fragments: {len(frag_df)}")


if __name__ == "__main__":
    main()
