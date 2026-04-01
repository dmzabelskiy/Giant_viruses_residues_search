#!/usr/bin/env python3
from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

try:
    from _pipeline import (
        RELAXED_CONSENSUS_CATEGORY,
        STRICT_CONSENSUS_CATEGORY,
        project_root,
    )
except Exception:
    import os
    DEFAULT_PROJECT_ROOT = Path.home() / "projects" / "gv_project"
    STRICT_CONSENSUS_CATEGORY = "rankpct_high_confidence"
    RELAXED_CONSENSUS_CATEGORY = "rankpct_moderate_confidence"
    def project_root() -> Path:
        return Path(os.environ.get("GV_PROJECT_ROOT", DEFAULT_PROJECT_ROOT))


def slugify_family_name(family: str) -> str:
    s = str(family).strip()
    s = re.sub(r"[^A-Za-z0-9]+", "_", s)
    s = re.sub(r"_+", "_", s).strip("_")
    return s.lower()


def family_stem(family: str, results_subdir: Optional[str] = None) -> str:
    return results_subdir or slugify_family_name(family)


def results_dir_for_family(family: str, results_subdir: Optional[str] = None) -> Path:
    return project_root() / "data" / "results" / family_stem(family, results_subdir)


def intermediate_dir_for_family(family: str, results_subdir: Optional[str] = None) -> Path:
    return results_dir_for_family(family, results_subdir) / "intermediate"


@dataclass(frozen=True)
class FamilyPaths:
    family: str
    stem: str
    results: Path
    intermediate: Path

    @property
    def main_ungapped_fasta(self) -> Path:
        return self.results / f"{self.stem}_main_ungapped.faa"

    @property
    def main_alignment_fasta(self) -> Path:
        return self.results / f"{self.stem}_main_mafft.trimmed.faa"

    @property
    def residue_calls_tsv(self) -> Path:
        return self.results / f"{self.stem}_viral_vs_cell_residue_calls.tsv"

    @property
    def residue_shortlist_tsv(self) -> Path:
        return self.results / f"{self.stem}_viral_vs_cell_residue_shortlist.tsv"

    @property
    def viral_full_length_index_tsv(self) -> Path:
        return self.results / f"viral_{self.stem}_full_length.tsv"

    @property
    def representative_ranking_tsv(self) -> Path:
        return self.results / f"viral_{self.stem}_representative_ranking.tsv"

    @property
    def representative_fasta(self) -> Path:
        return self.results / f"viral_{self.stem}_representative_for_colabfold.faa"

    @property
    def viral_raw_fasta(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_raw.faa"

    @property
    def viral_index_tsv(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_index.tsv"

    @property
    def viral_missing_tsv(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_missing.tsv"

    @property
    def viral_extraction_summary_tsv(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_extraction_summary.tsv"

    @property
    def alignment_membership_tsv(self) -> Path:
        return self.intermediate / f"{self.stem}_alignment_membership.tsv"

    @property
    def cellular_from_alignment_fasta(self) -> Path:
        return self.intermediate / f"cellular_{self.stem}_from_alignment.faa"

    @property
    def viral_from_alignment_fasta(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_from_alignment.faa"

    @property
    def alignment_unmatched_fasta(self) -> Path:
        return self.intermediate / f"{self.stem}_alignment_unmatched.faa"

    @property
    def alignment_summary_tsv(self) -> Path:
        return self.intermediate / f"{self.stem}_alignment_summary.tsv"

    @property
    def viral_full_length_fasta(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_full_length.faa"

    @property
    def viral_fragments_fasta(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_fragments.faa"

    @property
    def viral_fragments_tsv(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_fragments.tsv"

    @property
    def viral_length_filter_summary_tsv(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_length_filter_summary.tsv"

    @property
    def cellular_ungapped_fasta(self) -> Path:
        return self.intermediate / f"cellular_{self.stem}_ungapped.faa"

    @property
    def residue_calls_summary_tsv(self) -> Path:
        return self.intermediate / f"{self.stem}_viral_vs_cell_residue_calls_summary.tsv"

    @property
    def representative_excluded_tsv(self) -> Path:
        return self.intermediate / f"viral_{self.stem}_sequences_excluded_from_representative_due_to_X.tsv"

    @property
    def legacy_alignment_fasta(self) -> Path:
        return self.results / "aln" / f"{self.stem}_main_mafft.trimmed.faa"

    def ensure_dirs(self) -> None:
        self.results.mkdir(parents=True, exist_ok=True)
        self.intermediate.mkdir(parents=True, exist_ok=True)


def family_paths(family: str, results_subdir: Optional[str] = None) -> FamilyPaths:
    stem = family_stem(family, results_subdir)
    results = results_dir_for_family(family, results_subdir)
    return FamilyPaths(
        family=family,
        stem=stem,
        results=results,
        intermediate=intermediate_dir_for_family(family, results_subdir),
    )
