#!/usr/bin/env python3

from pathlib import Path
import pandas as pd
import re
from _pipeline import project_root

PROJECT_ROOT = project_root()
RAW = PROJECT_ROOT / "initial_raw"
INTERMEDIATE = PROJECT_ROOT / "data" / "intermediate"

INTERMEDIATE.mkdir(parents=True, exist_ok=True)


def normalize_mag_id(s: str) -> str:
    if pd.isna(s):
        return ""
    s = str(s).strip()
    s = re.sub(r"\.faa$", "", s)
    s = re.sub(r"\.fa\.dc$", "", s)
    s = re.sub(r"\.dc\.fa$", "", s)
    s = re.sub(r"\.fa$", "", s)
    return s


def pick_sheet(xls: pd.ExcelFile, candidates: list[str]) -> str:
    sheet_names = xls.sheet_names

    for cand in candidates:
        for s in sheet_names:
            if cand.lower() == s.lower():
                return s

    for cand in candidates:
        for s in sheet_names:
            if cand.lower() in s.lower():
                return s

    raise ValueError(f"Could not find a matching sheet among: {sheet_names}")


def main():
    s1_path = RAW / "Datasets" / "Dataset_S1.xlsx"
    s2_path = RAW / "Datasets" / "Dataset_S2.xlsx"
    prot_dir = RAW / "final_bins_prot"

    print(f"Reading workbook: {s1_path}")
    xls1 = pd.ExcelFile(s1_path)
    print("Dataset_S1 sheets:", xls1.sheet_names)

    print(f"\nReading workbook: {s2_path}")
    xls2 = pd.ExcelFile(s2_path)
    print("Dataset_S2 sheets:", xls2.sheet_names)

    # -----------------------------
    # Dataset S1: MAG metadata
    # -----------------------------
    s1_sheet = pick_sheet(xls1,["Final_MAGs_Statistics", "Final_MAG_stats", "MAG_stats", "Dataset 1", "Sheet1"])
    
    print(f"\nUsing Dataset_S1 sheet: {s1_sheet}")

    mag = pd.read_excel(s1_path, sheet_name=s1_sheet)
    mag.columns = [str(c).strip() for c in mag.columns]

    print("\nColumns in Dataset_S1 selected sheet:")
    print(list(mag.columns))

    if "MAG_ID" not in mag.columns:
        raise ValueError("Expected 'MAG_ID' column in Dataset_S1 selected sheet.")

    mag["mag_id_original"] = mag["MAG_ID"].astype(str)
    mag["mag_id_short"] = mag["mag_id_original"].map(normalize_mag_id)

    fasta_rows = []
    for fp in sorted(prot_dir.glob("*.faa")):
        fasta_rows.append({
            "fasta_file": fp.name,
            "fasta_path": str(fp.resolve()),
            "mag_id_short": normalize_mag_id(fp.name),
        })
    fasta_df = pd.DataFrame(fasta_rows)

    mag = pd.merge(mag, fasta_df, on="mag_id_short", how="left")

    # -----------------------------
    # Dataset S2: inspect sheets, but do not force OG parsing
    # -----------------------------
    s2_sheet = pick_sheet(
        xls2,
        ["Orthologous_Group_Annotation", "Orthologous Group Annotation", "Sheet1"]
    )
    print(f"\nUsing Dataset_S2 sheet: {s2_sheet}")

    og = pd.read_excel(s2_path, sheet_name=s2_sheet)
    og.columns = [str(c).strip() for c in og.columns]

    print("\nColumns in Dataset_S2 selected sheet:")
    print(list(og.columns))

    # Keep the raw sheet as-is. We only inspect/export it for now.
    # We already know the CSV export was not the per-OG table we wanted.
    # The workbook version is still worth preserving for later inspection.

    mag_out = INTERMEDIATE / "mag_metadata.tsv"
    fasta_out = INTERMEDIATE / "available_protein_fastas.tsv"
    s2_raw_out = INTERMEDIATE / "dataset_s2_selected_sheet.tsv"

    mag.to_csv(mag_out, sep="\t", index=False)
    fasta_df.to_csv(fasta_out, sep="\t", index=False)
    og.to_csv(s2_raw_out, sep="\t", index=False)

    print("\nWrote:")
    print(f"  {mag_out}")
    print(f"  {fasta_out}")
    print(f"  {s2_raw_out}")

    print("\nSummary:")
    print(f"  MAG metadata rows: {len(mag)}")
    print(f"  Protein FASTA files: {len(fasta_df)}")
    print(f"  Dataset_S2 selected sheet rows: {len(og)}")


if __name__ == "__main__":
    main()
