#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd


def standardize_group_label(raw: str) -> str:
    s = str(raw).strip()
    if not s:
        return "unknown"
    s_low = s.lower()
    if "this_study" in s_low or "this study" in s_low or "viral" in s_low:
        return "viral"
    return "cellular"


REQUIRED_COLUMNS = {"gene_family", "domain_group", "raw_label"}
OPTIONAL_COLUMNS = {"seq_core", "mag_id_short", "clade", "has_explicit_viral_suffix"}


def load_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Input table not found: {path}")

    df = pd.read_csv(path, sep="\t")
    missing = REQUIRED_COLUMNS - set(df.columns)
    if missing:
        raise ValueError(
            f"Input table is missing required columns: {sorted(missing)}\n"
            f"Found columns: {list(df.columns)}"
        )

    df = df.copy()
    df["gene_family"] = df["gene_family"].astype(str).str.strip()
    df["domain_group"] = df["domain_group"].astype(str).str.strip()
    df["raw_label"] = df["raw_label"].astype(str).str.strip()
    df["supergroup"] = df["domain_group"].map(standardize_group_label)

    # Prefer seq_core as the unique sequence identifier if available,
    # otherwise fall back to raw_label.
    if "seq_core" in df.columns:
        seq_id = df["seq_core"].fillna("").astype(str).str.strip()
        df["unique_seq_id"] = seq_id.where(seq_id.ne(""), df["raw_label"])
    else:
        df["unique_seq_id"] = df["raw_label"]

    if "mag_id_short" not in df.columns:
        df["mag_id_short"] = ""
    else:
        df["mag_id_short"] = df["mag_id_short"].fillna("").astype(str).str.strip()

    if "clade" not in df.columns:
        df["clade"] = ""
    else:
        df["clade"] = df["clade"].fillna("").astype(str).str.strip()

    return df


def family_summary(df: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    for family, sub in df.groupby("gene_family", dropna=False):
        viral = sub[sub["supergroup"] == "viral"].copy()
        cell = sub[sub["supergroup"] == "cellular"].copy()

        viral_unique = viral["unique_seq_id"].nunique()
        cell_unique = cell["unique_seq_id"].nunique()
        total_unique = sub["unique_seq_id"].nunique()

        rows.append(
            {
                "gene_family": family,
                "total_rows": len(sub),
                "total_unique_sequences": total_unique,
                "viral_rows": len(viral),
                "viral_unique_sequences": viral_unique,
                "viral_unique_mags": viral.loc[viral["mag_id_short"].ne(""), "mag_id_short"].nunique(),
                "viral_unique_clades": viral.loc[viral["clade"].ne(""), "clade"].nunique(),
                "cellular_rows": len(cell),
                "cellular_unique_sequences": cell_unique,
                "viral_fraction_of_unique_sequences": (viral_unique / total_unique) if total_unique else 0.0,
                "cellular_fraction_of_unique_sequences": (cell_unique / total_unique) if total_unique else 0.0,
            }
        )

    out = pd.DataFrame(rows)
    out = out.sort_values(
        ["viral_unique_sequences", "cellular_unique_sequences", "total_unique_sequences", "gene_family"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    return out


def domain_group_summary(df: pd.DataFrame) -> pd.DataFrame:
    out = (
        df.groupby(["gene_family", "domain_group", "supergroup"], dropna=False)
          .agg(
              rows=("raw_label", "size"),
              unique_sequences=("unique_seq_id", "nunique"),
              unique_mags=("mag_id_short", lambda s: s[s.ne("")].nunique()),
              unique_clades=("clade", lambda s: s[s.ne("")].nunique()),
          )
          .reset_index()
          .sort_values(["gene_family", "supergroup", "unique_sequences", "domain_group"], ascending=[True, True, False, True])
          .reset_index(drop=True)
    )
    return out


def candidate_families(summary: pd.DataFrame, min_viral: int, min_cellular: int) -> pd.DataFrame:
    cand = summary[
        (summary["viral_unique_sequences"] >= min_viral)
        & (summary["cellular_unique_sequences"] >= min_cellular)
    ].copy()

    cand["analysis_priority_score"] = (
        cand["viral_unique_sequences"].clip(upper=100)
        + 0.35 * cand["cellular_unique_sequences"].clip(upper=300)
        + 2.0 * cand["viral_unique_clades"]
    )

    cand = cand.sort_values(
        ["analysis_priority_score", "viral_unique_sequences", "cellular_unique_sequences", "gene_family"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)
    return cand


def main() -> None:
    ap = argparse.ArgumentParser(
        description="Summarize how many curated sequences are available for each gene family."
    )
    ap.add_argument(
        "--input",
        default="data/intermediate/curated_gene_sequences.tsv",
        help="TSV produced by 02_index_curated_gene_trees.py",
    )
    ap.add_argument(
        "--outdir",
        default="data/results/family_inventory",
        help="Directory for output TSV files",
    )
    ap.add_argument(
        "--min-viral",
        type=int,
        default=10,
        help="Minimum number of unique viral sequences to keep as an analysis candidate",
    )
    ap.add_argument(
        "--min-cellular",
        type=int,
        default=30,
        help="Minimum number of unique cellular sequences to keep as an analysis candidate",
    )
    args = ap.parse_args()

    input_path = Path(args.input).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_table(input_path)

    family = family_summary(df)
    domain = domain_group_summary(df)
    candidates = candidate_families(family, min_viral=args.min_viral, min_cellular=args.min_cellular)

    family_out = outdir / "gene_family_counts.tsv"
    domain_out = outdir / "gene_family_domain_group_counts.tsv"
    candidate_out = outdir / "gene_families_passing_minimums.tsv"

    family.to_csv(family_out, sep="\t", index=False)
    domain.to_csv(domain_out, sep="\t", index=False)
    candidates.to_csv(candidate_out, sep="\t", index=False)

    print(f"Input rows: {len(df)}")
    print(f"Unique gene families: {df['gene_family'].nunique()}")
    print()
    print(f"Wrote: {family_out}")
    print(f"Wrote: {domain_out}")
    print(f"Wrote: {candidate_out}")

    print("\nTop 20 families by viral unique-sequence count:")
    show_cols = [
        "gene_family",
        "viral_unique_sequences",
        "cellular_unique_sequences",
        "total_unique_sequences",
        "viral_unique_mags",
        "viral_unique_clades",
    ]
    print(family[show_cols].head(20).to_string(index=False))

    print(
        f"\nFamilies passing thresholds (viral >= {args.min_viral}, cellular >= {args.min_cellular}): {len(candidates)}"
    )
    if not candidates.empty:
        print(candidates[show_cols + ["analysis_priority_score"]].head(20).to_string(index=False))


if __name__ == "__main__":
    main()
