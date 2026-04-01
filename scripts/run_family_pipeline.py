#!/usr/bin/env python3
from __future__ import annotations

import argparse
import shlex
import shutil
import subprocess
import sys
from pathlib import Path

from _pipeline import project_root
from _family_pipeline import family_paths

PROJECT_ROOT = project_root()
SCRIPTS_DIR = PROJECT_ROOT / 'scripts'


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description=(
            'Run the generic family pipeline (steps 03-08) for any gene family. '
            'This wrapper orchestrates extraction, alignment-membership parsing, '
            'automatic length filtering, dataset preparation, MAFFT alignment, '
            'residue calling, and representative-sequence selection.'
        )
    )
    ap.add_argument('--family', required=True, help='Exact gene_family value from curated tables.')
    ap.add_argument('--alignment-source', required=True, help='Source family alignment FASTA used by step 04.')
    ap.add_argument('--results-subdir', default=None, help='Optional custom results directory name.')
    ap.add_argument('--skip-mafft', action='store_true', help='Skip MAFFT and reuse an existing aligned FASTA.')
    ap.add_argument('--alignment-final', default=None, help='Aligned FASTA to use for step 07. Defaults to results/<family>/<slug>_main_mafft.trimmed.faa')
    ap.add_argument('--mafft-bin', default='mafft', help='MAFFT executable name or path.')
    ap.add_argument('--mafft-mode', choices=['linsi', 'auto'], default='linsi', help='Alignment mode. linsi = --localpair --maxiterate 1000; auto = --auto.')
    ap.add_argument('--threads', type=int, default=None, help='Optional MAFFT --thread value.')
    ap.add_argument('--min-length', type=int, default=None, help='Manual override for viral full-length threshold.')
    ap.add_argument('--quantile', type=float, default=0.10, help='Auto-threshold lower cellular quantile.')
    ap.add_argument('--mad-scale', type=float, default=3.0, help='Auto-threshold robust MAD scale.')
    ap.add_argument('--stop-after', choices=['03', '04', '05', '06', 'mafft', '07', '08'], default=None, help='Stop after a specific stage.')
    ap.add_argument('--dry-run', action='store_true', help='Print commands without executing them.')
    return ap


def run_cmd(cmd: list[str], *, cwd: Path, dry_run: bool) -> None:
    printable = ' '.join(shlex.quote(x) for x in cmd)
    print(f'\n$ {printable}')
    if dry_run:
        return
    subprocess.run(cmd, cwd=str(cwd), check=True)


def ensure_exists(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f'{label} not found: {path}')


def script_path(name: str) -> Path:
    path = SCRIPTS_DIR / name
    ensure_exists(path, 'Required script')
    return path


def main() -> None:
    args = build_parser().parse_args()

    family = args.family
    paths = family_paths(family, args.results_subdir)
    paths.ensure_dirs()
    alignment_final = Path(args.alignment_final).expanduser().resolve() if args.alignment_final else paths.main_alignment_fasta
    canonical_alignment_final = paths.main_alignment_fasta
    alignment_source = Path(args.alignment_source).expanduser().resolve()

    ensure_exists(alignment_source, 'Source family alignment')
    if args.skip_mafft and not alignment_final.exists() and not args.dry_run:
        raise FileNotFoundError(f'--skip-mafft was set but aligned FASTA does not exist: {alignment_final}')
    if not args.skip_mafft and shutil.which(args.mafft_bin) is None and not args.dry_run:
        raise FileNotFoundError(f'Could not find MAFFT executable: {args.mafft_bin}')

    print('Family pipeline wrapper')
    print(f'  project root: {PROJECT_ROOT}')
    print(f'  family: {family}')
    print(f'  results dir: {paths.results}')
    print(f'  intermediate dir: {paths.intermediate}')
    print(f'  source alignment: {alignment_source}')
    print(f'  canonical final aligned FASTA: {canonical_alignment_final}')
    if alignment_final != canonical_alignment_final:
        print(f'  alignment path override: {alignment_final}')
    print(f'  MAFFT mode: {args.mafft_mode}')
    if args.min_length is None:
        print(f'  length threshold: auto (quantile={args.quantile}, mad_scale={args.mad_scale})')
    else:
        print(f'  length threshold: manual ({args.min_length} aa)')

    common = ['--family', family]
    if args.results_subdir:
        common += ['--results-subdir', args.results_subdir]

    # Step 03
    run_cmd([sys.executable, str(script_path('03_extract_family_sequences.py')), *common], cwd=PROJECT_ROOT, dry_run=args.dry_run)
    if args.stop_after == '03':
        return

    # Step 04
    run_cmd([
        sys.executable,
        str(script_path('04_extract_family_from_alignment.py')),
        *common,
        '--alignment',
        str(alignment_source),
    ], cwd=PROJECT_ROOT, dry_run=args.dry_run)
    if args.stop_after == '04':
        return

    # Step 05
    cmd05 = [
        sys.executable,
        str(script_path('05_filter_family_lengths.py')),
        *common,
        '--quantile', str(args.quantile),
        '--mad-scale', str(args.mad_scale),
    ]
    if args.min_length is not None:
        cmd05 += ['--min-length', str(args.min_length)]
    run_cmd(cmd05, cwd=PROJECT_ROOT, dry_run=args.dry_run)
    if args.stop_after == '05':
        return

    # Step 06
    run_cmd([sys.executable, str(script_path('06_prepare_family_main_dataset.py')), *common], cwd=PROJECT_ROOT, dry_run=args.dry_run)
    if args.stop_after == '06':
        return

    # MAFFT
    input_fasta = paths.main_ungapped_fasta
    ensure_exists(input_fasta, 'Combined ungapped FASTA') if not args.dry_run else None
    alignment_final.parent.mkdir(parents=True, exist_ok=True)
    canonical_alignment_final.parent.mkdir(parents=True, exist_ok=True)

    if not args.skip_mafft:
        mafft_cmd = [args.mafft_bin]
        if args.mafft_mode == 'linsi':
            mafft_cmd += ['--localpair', '--maxiterate', '1000']
        else:
            mafft_cmd += ['--auto']
        if args.threads is not None:
            mafft_cmd += ['--thread', str(args.threads)]
        mafft_cmd += [str(input_fasta)]

        printable = ' '.join(shlex.quote(x) for x in mafft_cmd) + ' > ' + shlex.quote(str(alignment_final))
        print(f'\n$ {printable}')
        if not args.dry_run:
            with open(alignment_final, 'w') as fh:
                subprocess.run(mafft_cmd, cwd=str(PROJECT_ROOT), check=True, stdout=fh)
            if alignment_final != canonical_alignment_final:
                shutil.copy2(alignment_final, canonical_alignment_final)
    elif not args.dry_run and alignment_final != canonical_alignment_final:
        ensure_exists(alignment_final, 'Aligned FASTA override')
        shutil.copy2(alignment_final, canonical_alignment_final)
    if args.stop_after == 'mafft':
        return

    # Step 07
    cmd07 = [
        sys.executable,
        str(script_path('07_call_family_residues.py')),
        *common,
        '--alignment',
        str(canonical_alignment_final),
    ]
    run_cmd(cmd07, cwd=PROJECT_ROOT, dry_run=args.dry_run)
    if args.stop_after == '07':
        return

    # Step 08
    run_cmd([sys.executable, str(script_path('08_choose_representative_viral_family.py')), *common], cwd=PROJECT_ROOT, dry_run=args.dry_run)

    print('\nPipeline complete.')
    print('Key outputs:')
    print(f'  main ungapped FASTA: {paths.main_ungapped_fasta}')
    print(f'  final alignment: {paths.main_alignment_fasta}')
    print(f'  residue calls: {paths.residue_calls_tsv}')
    print(f'  shortlist: {paths.residue_shortlist_tsv}')
    print(f'  viral full-length index: {paths.viral_full_length_index_tsv}')
    print(f'  representative ranking: {paths.representative_ranking_tsv}')
    print(f'  representative FASTA: {paths.representative_fasta}')


if __name__ == '__main__':
    main()
