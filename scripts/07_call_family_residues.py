#!/usr/bin/env python3
from __future__ import annotations

import argparse
import math
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd
from Bio import SeqIO
from scipy.stats import fisher_exact

from _family_pipeline import (
    RELAXED_CONSENSUS_CATEGORY,
    STRICT_CONSENSUS_CATEGORY,
    family_paths,
)

MAX_GAP_FRACTION = 0.20
STRICT_MIN_VIRAL_FREQ = 0.75
RELAXED_MIN_VIRAL_FREQ = 0.60
STRICT_MAX_QVALUE = 0.05
RELAXED_MAX_QVALUE = 0.20
STRICT_MIN_RANK_PERCENTILE_MEAN = 0.92
RELAXED_MIN_RANK_PERCENTILE_MEAN = 0.85

GAP_CHARS = {'-', '.'}
MISSING_CHARS = GAP_CHARS | {'X', '?'}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description="Call viral-vs-cellular residue signatures for an arbitrary family.")
    ap.add_argument("--family", required=True)
    ap.add_argument("--results-subdir", default=None)
    ap.add_argument("--alignment", default=None, help="Override aligned FASTA path. Default: results/<family>/<family>_main_mafft.trimmed.faa")
    return ap


def load_id_set(path: Path) -> set[str]:
    if not path.exists():
        raise FileNotFoundError(f'Missing FASTA: {path}')
    return {rec.id for rec in SeqIO.parse(str(path), 'fasta')}


def clean_residues(chars: Iterable[str]) -> List[str]:
    return [c for c in chars if c not in MISSING_CHARS]


def gap_fraction(chars: Iterable[str]) -> float:
    chars = list(chars)
    if not chars:
        return 1.0
    bad = sum(1 for c in chars if c in MISSING_CHARS)
    return bad / len(chars)


def counts_dict(chars: Iterable[str]) -> Dict[str, int]:
    return dict(Counter(clean_residues(chars)))


def major_residue(chars: Iterable[str]) -> Tuple[str, float, int, int]:
    counts = Counter(clean_residues(chars))
    if not counts:
        return '', 0.0, 0, 0
    aa, n = counts.most_common(1)[0]
    total = sum(counts.values())
    return aa, n / total, n, total


def counts_string(counts: Dict[str, int]) -> str:
    if not counts:
        return ''
    ordered = sorted(counts.items(), key=lambda kv: (-kv[1], kv[0]))
    return ';'.join(f'{aa}:{n}' for aa, n in ordered)


def freqs_dict(counts: Dict[str, int]) -> Dict[str, float]:
    total = sum(counts.values())
    if total == 0:
        return {}
    return {aa: n / total for aa, n in counts.items()}


def sequence_harmony_like(viral_counts: Dict[str, int], cell_counts: Dict[str, int]) -> float:
    pv = freqs_dict(viral_counts)
    pc = freqs_dict(cell_counts)
    overlap = 0.0
    for aa in set(pv) | set(pc):
        overlap += min(pv.get(aa, 0.0), pc.get(aa, 0.0))
    return 1.0 - overlap


def top_enriched_residue(viral_counts: Dict[str, int], cell_counts: Dict[str, int]) -> Dict[str, float | str | int]:
    pv = freqs_dict(viral_counts)
    pc = freqs_dict(cell_counts)
    residues = sorted(set(pv) | set(pc))

    if not residues:
        return {
            'aa': '',
            'viral_n': 0,
            'viral_total': 0,
            'cell_n': 0,
            'cell_total': 0,
            'viral_freq': 0.0,
            'cell_freq': 0.0,
            'freq_diff': 0.0,
            'odds_ratio': math.nan,
            'pvalue': math.nan,
        }

    viral_total = sum(viral_counts.values())
    cell_total = sum(cell_counts.values())
    best = None

    for aa in residues:
        v_n = viral_counts.get(aa, 0)
        c_n = cell_counts.get(aa, 0)
        v_f = pv.get(aa, 0.0)
        c_f = pc.get(aa, 0.0)
        diff = v_f - c_f
        record = {
            'aa': aa,
            'viral_n': v_n,
            'viral_total': viral_total,
            'cell_n': c_n,
            'cell_total': cell_total,
            'viral_freq': v_f,
            'cell_freq': c_f,
            'freq_diff': diff,
        }
        if best is None:
            best = record
        else:
            best_key = (best['freq_diff'], best['viral_freq'], -ord(str(best['aa'])[0]))
            rec_key = (record['freq_diff'], record['viral_freq'], -ord(str(record['aa'])[0]))
            if rec_key > best_key:
                best = record

    table = [
        [best['viral_n'], best['viral_total'] - best['viral_n']],
        [best['cell_n'], best['cell_total'] - best['cell_n']],
    ]
    odds_ratio, pvalue = fisher_exact(table, alternative='two-sided')
    best['odds_ratio'] = odds_ratio
    best['pvalue'] = pvalue
    return best


def benjamini_hochberg(pvalues: List[float]) -> List[float]:
    indexed = [(i, p) for i, p in enumerate(pvalues)]
    valid = [(i, p) for i, p in indexed if pd.notna(p)]
    qvals = [math.nan] * len(pvalues)
    if not valid:
        return qvals

    valid_sorted = sorted(valid, key=lambda x: x[1])
    m = len(valid_sorted)
    adjusted = [0.0] * m
    prev = 1.0
    for rank in range(m, 0, -1):
        i, p = valid_sorted[rank - 1]
        q = min(prev, p * m / rank)
        adjusted[rank - 1] = q
        prev = q
    for (idx, _), q in zip(valid_sorted, adjusted):
        qvals[idx] = q
    return qvals


def assign_rank(series: pd.Series, ascending: bool) -> pd.Series:
    return series.rank(method='min', ascending=ascending)


def rank_to_percentile(rank: float, n: int) -> float:
    if pd.isna(rank) or n <= 0:
        return math.nan
    if n == 1:
        return 1.0
    return 1.0 - (float(rank) - 1.0) / float(n - 1)


def classify_consensus(row: pd.Series) -> str:
    if not row['passes_gap_filter']:
        return ''
    if row['top_viral_residue_freq'] <= row['top_cell_residue_freq']:
        return ''

    if (
        row['top_viral_residue_freq'] >= STRICT_MIN_VIRAL_FREQ
        and pd.notna(row['fisher_qvalue']) and row['fisher_qvalue'] <= STRICT_MAX_QVALUE
        and pd.notna(row['rank_percentile_mean']) and row['rank_percentile_mean'] >= STRICT_MIN_RANK_PERCENTILE_MEAN
    ):
        return STRICT_CONSENSUS_CATEGORY

    if (
        row['top_viral_residue_freq'] >= RELAXED_MIN_VIRAL_FREQ
        and pd.notna(row['fisher_qvalue']) and row['fisher_qvalue'] <= RELAXED_MAX_QVALUE
        and pd.notna(row['rank_percentile_mean']) and row['rank_percentile_mean'] >= RELAXED_MIN_RANK_PERCENTILE_MEAN
    ):
        return RELAXED_CONSENSUS_CATEGORY

    return ''


def main() -> None:
    args = build_parser().parse_args()
    paths = family_paths(args.family, args.results_subdir)
    paths.ensure_dirs()
    if args.alignment:
        align_path = Path(args.alignment).expanduser().resolve()
    elif paths.main_alignment_fasta.exists():
        align_path = paths.main_alignment_fasta
    else:
        align_path = paths.legacy_alignment_fasta
    viral_fasta = paths.viral_full_length_fasta
    cell_fasta = paths.cellular_ungapped_fasta

    out_all = paths.residue_calls_tsv
    out_summary = paths.residue_calls_summary_tsv
    out_shortlist = paths.residue_shortlist_tsv

    viral_ids = load_id_set(viral_fasta)
    cell_ids = load_id_set(cell_fasta)
    overlap = viral_ids & cell_ids
    if overlap:
        raise ValueError(f'IDs overlap between viral and cellular source FASTAs: {len(overlap)}')

    records = list(SeqIO.parse(str(align_path), 'fasta'))
    if not records:
        raise ValueError(f'No records found in {align_path}')

    aln_len = len(records[0].seq)
    for rec in records:
        if len(rec.seq) != aln_len:
            raise ValueError('Alignment is not rectangular.')

    seq_info = []
    unknown_ids = []
    for rec in records:
        if rec.id in viral_ids:
            group = 'viral'
        elif rec.id in cell_ids:
            group = 'cellular'
        else:
            group = 'unknown'
            unknown_ids.append(rec.id)
        seq_info.append({'seq_id': rec.id, 'description': rec.description, 'group': group, 'seq': str(rec.seq)})

    if unknown_ids:
        preview = ', '.join(unknown_ids[:10])
        raise ValueError(f'{len(unknown_ids)} aligned records were not found in viral/cellular source FASTAs. First few: {preview}')

    df = pd.DataFrame(seq_info)
    viral = df[df['group'] == 'viral'].copy()
    cell = df[df['group'] == 'cellular'].copy()

    print('Sequence counts:')
    print(f'  viral: {len(viral)}')
    print(f'  cellular: {len(cell)}')

    rows = []
    for pos in range(aln_len):
        viral_chars = [seq[pos] for seq in viral['seq']]
        cell_chars = [seq[pos] for seq in cell['seq']]

        viral_counts = counts_dict(viral_chars)
        cell_counts = counts_dict(cell_chars)
        viral_major, viral_major_freq, viral_major_n, viral_total = major_residue(viral_chars)
        cell_major, cell_major_freq, cell_major_n, cell_total = major_residue(cell_chars)
        top = top_enriched_residue(viral_counts, cell_counts)
        sh_like = sequence_harmony_like(viral_counts, cell_counts)
        v_gap = gap_fraction(viral_chars)
        c_gap = gap_fraction(cell_chars)

        rows.append({
            'alignment_pos_1based': pos + 1,
            'viral_consensus': viral_major,
            'viral_consensus_freq': viral_major_freq,
            'viral_consensus_support': f'{viral_major_n}/{viral_total}' if viral_total else '',
            'viral_non_gap_n': viral_total,
            'viral_gap_fraction': v_gap,
            'viral_counts': counts_string(viral_counts),
            'cell_consensus': cell_major,
            'cell_consensus_freq': cell_major_freq,
            'cell_consensus_support': f'{cell_major_n}/{cell_total}' if cell_total else '',
            'cell_non_gap_n': cell_total,
            'cell_gap_fraction': c_gap,
            'cell_counts': counts_string(cell_counts),
            'sequence_harmony_like': sh_like,
            'top_viral_residue': top['aa'],
            'top_viral_residue_support': f"{top['viral_n']}/{top['viral_total']}" if top['viral_total'] else '',
            'top_cell_residue_support': f"{top['cell_n']}/{top['cell_total']}" if top['cell_total'] else '',
            'top_viral_residue_freq': top['viral_freq'],
            'top_cell_residue_freq': top['cell_freq'],
            'top_freq_diff': top['freq_diff'],
            'odds_ratio': top['odds_ratio'],
            'fisher_pvalue': top['pvalue'],
            'passes_gap_filter': (v_gap <= MAX_GAP_FRACTION) and (c_gap <= MAX_GAP_FRACTION),
        })

    out = pd.DataFrame(rows)
    out['fisher_qvalue'] = benjamini_hochberg(out['fisher_pvalue'].tolist())

    eligible = (
        out['passes_gap_filter']
        & out['top_viral_residue'].astype(str).ne('')
        & out['top_freq_diff'].gt(0)
        & out['fisher_qvalue'].notna()
    )

    out['rank_q'] = math.nan
    out['rank_harmony'] = math.nan
    out['rank_freq_diff'] = math.nan
    out['rank_q_pct'] = math.nan
    out['rank_harmony_pct'] = math.nan
    out['rank_freq_diff_pct'] = math.nan
    out['rank_percentile_mean'] = math.nan
    out.loc[eligible, 'rank_q'] = assign_rank(out.loc[eligible, 'fisher_qvalue'], ascending=True)
    out.loc[eligible, 'rank_harmony'] = assign_rank(out.loc[eligible, 'sequence_harmony_like'], ascending=False)
    out.loc[eligible, 'rank_freq_diff'] = assign_rank(out.loc[eligible, 'top_freq_diff'], ascending=False)
    n_eligible = int(eligible.sum())
    out.loc[eligible, 'rank_q_pct'] = out.loc[eligible, 'rank_q'].map(lambda r: rank_to_percentile(r, n_eligible))
    out.loc[eligible, 'rank_harmony_pct'] = out.loc[eligible, 'rank_harmony'].map(lambda r: rank_to_percentile(r, n_eligible))
    out.loc[eligible, 'rank_freq_diff_pct'] = out.loc[eligible, 'rank_freq_diff'].map(lambda r: rank_to_percentile(r, n_eligible))
    out.loc[eligible, 'rank_percentile_mean'] = out.loc[
        eligible, ['rank_q_pct', 'rank_harmony_pct', 'rank_freq_diff_pct']
    ].mean(axis=1)

    # Retain legacy fixed-top-N columns as secondary diagnostics only.
    out['legacy_n_top15_axes'] = (
        out[['rank_q', 'rank_harmony', 'rank_freq_diff']].le(15).sum(axis=1).where(eligible, other=math.nan)
    )
    out['legacy_n_top30_axes'] = (
        out[['rank_q', 'rank_harmony', 'rank_freq_diff']].le(30).sum(axis=1).where(eligible, other=math.nan)
    )
    out['category'] = out.apply(classify_consensus, axis=1)
    out['ranking_score'] = out['rank_percentile_mean']
    out['site_weight_score'] = (
        out['rank_percentile_mean'].fillna(0.0)
        * out['top_viral_residue_freq'].clip(lower=0, upper=1).fillna(0.0)
    )

    out = out.sort_values(['category', 'ranking_score', 'alignment_pos_1based'], ascending=[True, False, True])
    out.to_csv(out_all, sep='\t', index=False)

    shortlist = out[out['category'] != ''].copy()
    shortlist = shortlist.sort_values(
        ['category', 'rank_percentile_mean', 'rank_q_pct', 'rank_harmony_pct', 'rank_freq_diff_pct', 'alignment_pos_1based'],
        ascending=[True, False, False, False, False, True],
    )
    shortlist.to_csv(out_shortlist, sep='\t', index=False)

    summary = pd.DataFrame({
        'metric': [
            'family', 'n_alignment_sequences', 'n_viral_sequences', 'n_cellular_sequences', 'n_positions',
            'n_passes_gap_filter', 'n_eligible_positions', f'n_{STRICT_CONSENSUS_CATEGORY}', f'n_{RELAXED_CONSENSUS_CATEGORY}'
        ],
        'value': [
            args.family, len(df), len(viral), len(cell), len(out), int(out['passes_gap_filter'].sum()), n_eligible,
            int((out['category'] == STRICT_CONSENSUS_CATEGORY).sum()),
            int((out['category'] == RELAXED_CONSENSUS_CATEGORY).sum()),
        ],
    })
    summary.to_csv(out_summary, sep='\t', index=False)

    print(f'\nWrote: {out_all}')
    print(f'Wrote: {out_shortlist}')
    print(f'Wrote: {out_summary}')
    print('\nCategory counts:')
    print(out['category'].replace('', 'uncategorized').value_counts().to_string())


if __name__ == '__main__':
    main()
