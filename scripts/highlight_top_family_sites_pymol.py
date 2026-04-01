#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path
from typing import Iterable

import pandas as pd
from Bio import SeqIO

GAPLIKE = set('-.')
MISSING = set('-.X?')


def compress_resi(nums: list[int]) -> str:
    if not nums:
        return ''
    nums = sorted(set(nums))
    chunks: list[str] = []
    start = prev = nums[0]
    for n in nums[1:]:
        if n == prev + 1:
            prev = n
            continue
        chunks.append(f'{start}-{prev}' if start != prev else str(start))
        start = prev = n
    chunks.append(f'{start}-{prev}' if start != prev else str(start))
    return '+'.join(chunks)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(
        description='Write a PyMOL script that highlights shortlisted family-specific sites after mapping alignment columns to model residue numbers.'
    )
    ap.add_argument('--shortlist', required=True, help='Shortlist TSV from 07_call_family_residues.py')
    ap.add_argument('--alignment', required=True, help='Final alignment FASTA used for residue calling')
    ap.add_argument('--representative-id', required=True, help='Sequence ID of the representative used for the structure model')
    ap.add_argument('--pdb', required=True, help='Representative structure model (PDB)')
    ap.add_argument('--init', default=None, help='Optional init.py style script to run inside PyMOL')
    ap.add_argument('--output', default=None, help='Output PyMOL .py path')
    ap.add_argument('--mapping-out', default=None, help='Optional mapped-site TSV output path')
    ap.add_argument('--top-n', type=int, default=15, help='Number of top sites to highlight after sorting')
    ap.add_argument('--category', default='auto', help='Category filter. Use auto to prefer strict if present, otherwise all non-empty categories.')
    ap.add_argument('--color', default='orange', help='PyMOL color name for highlighted sites when not splitting by confidence.')
    ap.add_argument('--split-by-confidence', action='store_true', help='Create separate PyMOL selections/colors for high- and moderate-confidence rows.')
    ap.add_argument('--high-color', default='tv_orange', help='PyMOL color for rankpct_high_confidence sites when --split-by-confidence is used.')
    ap.add_argument('--moderate-color', default='marine', help='PyMOL color for rankpct_moderate_confidence sites when --split-by-confidence is used.')
    ap.add_argument('--object-name', default='family_model', help='PyMOL object name')
    ap.add_argument('--selection-name', default='top_family_sites', help='PyMOL selection name')
    ap.add_argument('--chain', default=None, help='Optional PDB chain to use. Default: first chain with CA atoms.')
    ap.add_argument('--pdb-start-resi', type=int, default=None, help='Override PDB numbering by mapping ungapped position 1 to this residue number. Use only for simple continuous numbering.')
    return ap.parse_args()


def load_alignment_sequence(alignment_path: Path, representative_id: str) -> str:
    for rec in SeqIO.parse(str(alignment_path), 'fasta'):
        if rec.id == representative_id:
            return str(rec.seq)
    raise ValueError(f'Representative ID not found in alignment: {representative_id}')


def alignment_col_to_ungapped_index(aligned_seq: str) -> dict[int, int]:
    mapping: dict[int, int] = {}
    ungapped_idx = 0
    for i, aa in enumerate(aligned_seq, start=1):
        if aa not in GAPLIKE:
            ungapped_idx += 1
            mapping[i] = ungapped_idx
    return mapping


def parse_pdb_residue_order(pdb_path: Path, chain: str | None) -> tuple[str, list[int]]:
    residues: list[tuple[str, int, str]] = []
    seen: set[tuple[str, int, str]] = set()
    chosen_chain = chain

    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            atom_name = line[12:16].strip()
            if atom_name != 'CA':
                continue
            resname = line[17:20].strip()
            chain_id = line[21].strip() or '_'
            resseq_text = line[22:26].strip()
            icode = line[26].strip()
            if not resseq_text:
                continue
            try:
                resseq = int(resseq_text)
            except ValueError:
                continue
            if resname in {'HOH', 'WAT'}:
                continue
            if chosen_chain is None:
                chosen_chain = chain_id
            if chain_id != chosen_chain:
                continue
            key = (chain_id, resseq, icode)
            if key not in seen:
                seen.add(key)
                residues.append(key)

    if chosen_chain is None or not residues:
        raise ValueError(f'No CA residues found in PDB: {pdb_path}')

    if any(k[2] for k in residues):
        raise ValueError(
            'PDB contains insertion codes in the selected chain. Please renumber the model or use a simpler numbering scheme before generating PyMOL selections.'
        )

    return chosen_chain, [k[1] for k in residues]


def choose_sites(shortlist: pd.DataFrame, top_n: int, category: str) -> tuple[pd.DataFrame, str]:
    if 'alignment_pos_1based' not in shortlist.columns:
        raise ValueError('Shortlist TSV is missing alignment_pos_1based.')

    df = shortlist.copy()
    if 'category' in df.columns:
        df['category'] = df['category'].fillna('')
    else:
        df['category'] = ''

    if category == 'auto':
        cats = set(df['category']) - {''}
        if 'rankpct_high_confidence' in cats:
            chosen_category = 'rankpct_high_confidence'
            df = df[df['category'] == chosen_category].copy()
        else:
            chosen_category = 'all_nonempty'
            df = df[df['category'] != ''].copy()
    elif category.lower() == 'all':
        chosen_category = 'all_rows'
    elif category.lower() == 'all_nonempty':
        chosen_category = 'all_nonempty'
        df = df[df['category'] != ''].copy()
    else:
        chosen_category = category
        df = df[df['category'] == category].copy()

    if df.empty:
        raise ValueError(f'No rows remain after category filter: {chosen_category}')

    sort_cols: list[str] = []
    ascending: list[bool] = []
    for col in ['rank_percentile_mean', 'rank_q_pct', 'rank_harmony_pct', 'rank_freq_diff_pct', 'ranking_score', 'alignment_pos_1based']:
        if col in df.columns:
            sort_cols.append(col)
            ascending.append(col not in {'rank_percentile_mean', 'rank_q_pct', 'rank_harmony_pct', 'rank_freq_diff_pct', 'ranking_score'})
    if sort_cols:
        df = df.sort_values(sort_cols, ascending=ascending)

    return df.head(top_n).copy(), chosen_category


def build_pymol_script(
    pdb_path: Path,
    init_path: Path | None,
    object_name: str,
    selection_name: str,
    mapped_df: pd.DataFrame,
    color: str,
    split_by_confidence: bool,
    high_color: str,
    moderate_color: str,
    chain: str,
) -> str:
    resi_expr = compress_resi([int(x) for x in mapped_df['pdb_resi'].tolist()])
    lines: list[str] = [
        'from pymol import cmd',
        '',
        'cmd.reinitialize()',
    ]
    if init_path is not None:
        lines.append(f'cmd.do("run {init_path.as_posix()}")')
    lines.extend([
        f'cmd.load("{pdb_path.as_posix()}", "{object_name}")',
        f'cmd.hide("everything", "{object_name}")',
        f'cmd.show("cartoon", "{object_name}")',
        f'cmd.color("gray80", "{object_name}")',
        f'cmd.set("cartoon_transparency", 0.15, "{object_name}")',
    ])

    lines.extend([
        f'cmd.select("{selection_name}", "{object_name} and chain {chain} and resi {resi_expr}")',
        f'cmd.show("sticks", "{selection_name}")',
        f'cmd.show("spheres", "{selection_name} and name CA")',
        f'cmd.set("sphere_scale", 0.34, "{selection_name} and name CA")',
        f'cmd.set("stick_radius", 0.24, "{selection_name}")',
    ])

    if split_by_confidence and 'category' in mapped_df.columns:
        for category, suffix, cat_color in [
            ('rankpct_high_confidence', 'high_confidence', high_color),
            ('rankpct_moderate_confidence', 'moderate_confidence', moderate_color),
        ]:
            cat_rows = mapped_df[mapped_df['category'] == category].copy()
            if cat_rows.empty:
                continue
            cat_resi_expr = compress_resi([int(x) for x in cat_rows['pdb_resi'].tolist()])
            cat_selection = f'{selection_name}_{suffix}'
            lines.extend([
                f'cmd.select("{cat_selection}", "{object_name} and chain {chain} and resi {cat_resi_expr}")',
                f'cmd.show("sticks", "{cat_selection}")',
                f'cmd.show("spheres", "{cat_selection} and name CA")',
                f'cmd.color("{cat_color}", "{cat_selection}")',
            ])
    else:
        lines.append(f'cmd.color("{color}", "{selection_name}")')

    lines.extend([
        f'cmd.orient("{selection_name}")',
        f'cmd.zoom("{selection_name}", 10)',
        ''
    ])
    return '\n'.join(lines)


def main() -> None:
    args = parse_args()

    shortlist_path = Path(args.shortlist).expanduser().resolve()
    alignment_path = Path(args.alignment).expanduser().resolve()
    pdb_path = Path(args.pdb).expanduser().resolve()
    init_path = Path(args.init).expanduser().resolve() if args.init else None

    output_path = Path(args.output).expanduser().resolve() if args.output else pdb_path.with_name(pdb_path.stem + f'_top{args.top_n}_family_sites.py')
    mapping_out = Path(args.mapping_out).expanduser().resolve() if args.mapping_out else output_path.with_suffix('.mapped_sites.tsv')

    shortlist = pd.read_csv(shortlist_path, sep='\t')
    top, chosen_category = choose_sites(shortlist, args.top_n, args.category)

    aligned_seq = load_alignment_sequence(alignment_path, args.representative_id)
    ungapped_seq = ''.join(aa for aa in aligned_seq if aa not in GAPLIKE)
    col_to_seq = alignment_col_to_ungapped_index(aligned_seq)

    if args.pdb_start_resi is not None:
        chain_id = args.chain or 'A'
        pdb_resnums = [args.pdb_start_resi + i for i in range(len(ungapped_seq))]
    else:
        chain_id, pdb_resnums = parse_pdb_residue_order(pdb_path, args.chain)
        if len(pdb_resnums) != len(ungapped_seq):
            raise ValueError(
                'Representative ungapped alignment length does not match the number of CA residues in the selected PDB chain. '
                f'Ungapped alignment length={len(ungapped_seq)}, PDB CA residues={len(pdb_resnums)}, chain={chain_id}. '
                'Use the exact representative model or provide --pdb-start-resi only when numbering is simple and continuous.'
            )

    mapped_rows: list[dict[str, object]] = []
    skipped_rows: list[dict[str, object]] = []
    for _, row in top.iterrows():
        aln_pos = int(row['alignment_pos_1based'])
        if aln_pos < 1 or aln_pos > len(aligned_seq):
            skipped_rows.append({'alignment_pos_1based': aln_pos, 'reason': 'outside_alignment'})
            continue

        rep_aa = aligned_seq[aln_pos - 1]
        if rep_aa in MISSING:
            skipped_rows.append({'alignment_pos_1based': aln_pos, 'reason': f'representative_has_{rep_aa}'})
            continue

        seq_idx = col_to_seq.get(aln_pos)
        if seq_idx is None:
            skipped_rows.append({'alignment_pos_1based': aln_pos, 'reason': 'failed_alignment_to_sequence_mapping'})
            continue

        pdb_resi = pdb_resnums[seq_idx - 1]
        out_row = row.to_dict()
        out_row['representative_aligned_aa'] = rep_aa
        out_row['representative_seq_pos_1based'] = seq_idx
        out_row['pdb_chain'] = chain_id
        out_row['pdb_resi'] = pdb_resi
        mapped_rows.append(out_row)

    if not mapped_rows:
        raise ValueError('No shortlisted sites could be mapped from alignment columns to PDB residue numbers.')

    mapped_df = pd.DataFrame(mapped_rows)
    mapped_df.to_csv(mapping_out, sep='\t', index=False)

    resi = [int(x) for x in mapped_df['pdb_resi'].tolist()]
    output_path.write_text(
        build_pymol_script(
            pdb_path,
            init_path,
            args.object_name,
            args.selection_name,
            mapped_df,
            args.color,
            args.split_by_confidence,
            args.high_color,
            args.moderate_color,
            chain_id,
        )
    )

    print(f'Shortlist: {shortlist_path}')
    print(f'Alignment: {alignment_path}')
    print(f'Representative: {args.representative_id}')
    print(f'PDB: {pdb_path}')
    print(f'Chosen category filter: {chosen_category}')
    print(f'Selected rows before mapping: {len(top)}')
    print(f'Mapped rows: {len(mapped_df)}')
    print(f'Skipped rows: {len(skipped_rows)}')
    if skipped_rows:
        print('Skipped alignment positions: ' + ', '.join(str(x['alignment_pos_1based']) for x in skipped_rows))
    print(f'PDB chain used: {chain_id}')
    print(f'Wrote mapped-site table: {mapping_out}')
    print(f'Wrote PyMOL script: {output_path}')
    print('Highlighted PDB residues: ' + ', '.join(str(x) for x in sorted(set(resi))))
    print('\nOpen with:')
    print(f'pymol -r {output_path.as_posix()}')


if __name__ == '__main__':
    main()
