#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
import re
import pandas as pd
from _pipeline import project_root

PROJECT_ROOT = project_root()
RAW = PROJECT_ROOT / "initial_raw"
INTERMEDIATE = PROJECT_ROOT / "data" / "intermediate"
MARKER_DIR = RAW / "marker_gene_phylogenies" / "Phylogenetic_tree_alignments_and_tree_files"
XLS_PATH = MARKER_DIR / "Domain_level_classification_of_tree_sequences.xls"

CLADES = ("MM", "EP", "LP", "IR", "MR", "PT", "PX", "AF")
VIRAL_SUFFIX_RE = re.compile(r"\.\.\.((?:MM|EP|LP|IR|MR|PT|PX|AF)\d+)\s*$")
MAG_PREFIX_RE = re.compile(r"((?:ERX|SRX)\d+\.\d+)")


def normalize_mag_id_from_label(label: str) -> str:
    m = MAG_PREFIX_RE.search(str(label))
    return m.group(1) if m else ""


def extract_clade_from_label(label: str) -> str:
    """Only trust clades when they appear as the explicit viral suffix `...MM01`."""
    s = str(label).strip()
    m = VIRAL_SUFFIX_RE.search(s)
    return m.group(1) if m else ""


def extract_seq_core(label: str) -> str:
    s = str(label).strip()
    m = re.match(r"^((?:ERX|SRX)\d+\.\d+\.fa\.dc\.\.[^.]+(?:_[^.]+)*)", s)
    if m:
        return m.group(1)
    return VIRAL_SUFFIX_RE.sub("", s)


def looks_like_sequence_label(value: str) -> bool:
    s = str(value).strip()
    if not s:
        return False
    lowered = s.lower()
    if lowered in {"nan", "none"}:
        return False
    # Keep broad coverage, but skip obvious placeholders/notes.
    if len(s) < 3:
        return False
    if lowered.startswith("unnamed:"):
        return False
    return True


def main() -> None:
    if not XLS_PATH.exists():
        raise FileNotFoundError(f"Could not find workbook: {XLS_PATH}")

    INTERMEDIATE.mkdir(parents=True, exist_ok=True)
    print(f"Using workbook: {XLS_PATH}")
    xls = pd.ExcelFile(XLS_PATH)
    print("Sheets:", xls.sheet_names)

    rows: list[dict[str, str]] = []
    skipped = 0
    for sheet in xls.sheet_names:
        df = pd.read_excel(XLS_PATH, sheet_name=sheet)
        if df.empty:
            continue
        df.columns = [str(c).strip() for c in df.columns]
        for col in df.columns:
            vals = df[col].dropna().astype(str)
            for raw in vals:
                v = raw.strip()
                if not looks_like_sequence_label(v):
                    skipped += 1
                    continue
                rows.append(
                    {
                        "gene_family": sheet,
                        "domain_group": col,
                        "raw_label": v,
                        "seq_core": extract_seq_core(v),
                        "mag_id_short": normalize_mag_id_from_label(v),
                        "clade": extract_clade_from_label(v),
                        "has_explicit_viral_suffix": bool(VIRAL_SUFFIX_RE.search(v)),
                    }
                )

    out = pd.DataFrame(rows)
    out = out.drop_duplicates(subset=["gene_family", "domain_group", "raw_label"]).reset_index(drop=True)

    viral_mask = out["domain_group"].str.contains("This_study", case=False, na=False)
    suspicious = out[(out["clade"] != "") & (~viral_mask)]

    all_out = INTERMEDIATE / "curated_gene_sequences.tsv"
    viral_out = INTERMEDIATE / "curated_viral_sequences.tsv"
    cellular_out = INTERMEDIATE / "curated_cellular_sequences.tsv"
    suspicious_out = INTERMEDIATE / "curated_gene_sequences_suspicious_clades.tsv"

    out.to_csv(all_out, sep="\t", index=False)
    out.loc[viral_mask].to_csv(viral_out, sep="\t", index=False)
    out.loc[~viral_mask].to_csv(cellular_out, sep="\t", index=False)
    suspicious.to_csv(suspicious_out, sep="\t", index=False)

    print("\nWrote:")
    print(f"  {all_out}")
    print(f"  {viral_out}")
    print(f"  {cellular_out}")
    print(f"  {suspicious_out}")

    print("\nCounts by family/domain:")
    summary = out.groupby(["gene_family", "domain_group"]).size().reset_index(name="n")
    print(summary.to_string(index=False))
    print(f"\nIndexed rows: {len(out)}")
    print(f"Skipped obvious non-label cells: {skipped}")
    print(f"Non-viral rows carrying explicit clade suffix: {len(suspicious)}")


if __name__ == "__main__":
    main()
