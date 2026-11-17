#!/usr/bin/env python
"""Compare package outputs with reference outputs (for example the notebooks').

Usage
-----
1. Copy the outputs of the notebooks to a reference folder (the package writes
   to the same places and would overwrite them)::

       mkdir reference
       cp <root>/Tables_*/Integrated_*/*_longitudinal.csv reference/
       cp <root>/Integrated_ALL/*.csv reference/

2. Run the package (``inbrain-dementia preprocess <root>`` and
   ``inbrain-dementia integrate <root>``; cohort tools with the same settings).

3. Compare::

       python tools/compare_outputs.py <root> reference/ [--cohorts-dir <folder>]

Every reference file with a known name is compared with the package file of
the same name: same rows (by key), same columns (order ignored), same values
(numbers compared after rounding to 1e-9, so ``1`` equals ``1.0``). Wide files
are compared in blocks of columns to bound memory.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from inbrain_dementia_data_analyser.config import load_config  # noqa: E402

BLOCK = 1500


def _read(path, key, cols):
    parts = pd.read_csv(path, usecols=[key] + cols, dtype=str, keep_default_na=False, chunksize=5000)
    return pd.concat(parts, ignore_index=True)


def compare(new: Path, ref: Path, key: str | None) -> list[str]:
    """Return a list of problems (empty when the files agree)."""
    hn = pd.read_csv(new, nrows=0).columns.tolist()
    hr = pd.read_csv(ref, nrows=0).columns.tolist()
    problems = []
    if set(hn) != set(hr):
        problems.append(f"columns only in new: {sorted(set(hn) - set(hr))[:10]}; "
                        f"only in reference: {sorted(set(hr) - set(hn))[:10]}")
    common = [c for c in hr if c in set(hn) and c != key]
    if key is None:
        key = "__row"
    for a in range(0, len(common), BLOCK):
        cols = common[a:a + BLOCK]
        if key == "__row":
            x = pd.read_csv(new, usecols=cols, dtype=str, keep_default_na=False)
            y = pd.read_csv(ref, usecols=cols, dtype=str, keep_default_na=False)
            if len(x) != len(y):
                return problems + [f"row count {len(x)} vs {len(y)}"]
        else:
            x = _read(new, key, cols).drop_duplicates(key).set_index(key)
            y = _read(ref, key, cols).drop_duplicates(key).set_index(key)
            if set(x.index) != set(y.index):
                problems.append(f"keys only in new: {len(set(x.index) - set(y.index))}, "
                                f"only in reference: {len(set(y.index) - set(x.index))}")
            common_idx = y.index.intersection(x.index)
            x, y = x.loc[common_idx], y.loc[common_idx]
        for c in cols:
            d = x[c].values != y[c].values
            if not d.any():
                continue
            nx = pd.to_numeric(pd.Series(x[c].values[d]), errors="coerce")
            ny = pd.to_numeric(pd.Series(y[c].values[d]), errors="coerce")
            real = ~((nx - ny).abs() < 1e-9)
            if real.any():
                ex = list(zip(pd.Series(x[c].values[d])[real].head(2), pd.Series(y[c].values[d])[real].head(2)))
                problems.append(f"{c}: {int(real.sum())} values differ, e.g. {ex}")
    return problems


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("root", help="root folder the package wrote to")
    ap.add_argument("reference", help="folder with the reference CSV files")
    ap.add_argument("--cohorts-dir", default=None,
                    help="folder with the package's cohort-tool outputs (default: <root>/Integrated_ALL)")
    ap.add_argument("--config-dir", default=None)
    args = ap.parse_args(argv)
    cfg = load_config(args.config_dir)
    root, ref = Path(args.root), Path(args.reference)
    integ = root / cfg.integration.get("output_folder", "Integrated_ALL")
    cohorts_dir = Path(args.cohorts_dir) if args.cohorts_dir else integ
    targets = {ds.output_file: (ds.output_path(root), "ID") for ds in cfg.datasets.values()}
    targets["MASTER_all_cohorts.csv"] = (integ / "MASTER_all_cohorts.csv", "uid")
    targets["MASTER_column_coverage.csv"] = (integ / "MASTER_column_coverage.csv", "Unnamed: 0")
    targets["MASTER_expanded_240m.csv"] = (integ / "MASTER_expanded_240m.csv", "uid")
    n_ok = n_bad = 0
    for rf in sorted(ref.glob("*.csv")):
        if rf.name in targets:
            new, key = targets[rf.name]
        else:
            new, key = cohorts_dir / rf.name, None     # cohort outputs: compared row by row
        if not new.exists():
            print(f"[missing] {rf.name}: no package output at {new}")
            n_bad += 1
            continue
        probs = compare(new, rf, key)
        if probs:
            n_bad += 1
            print(f"[DIFF] {rf.name}")
            for p in probs[:12]:
                print(f"    {p}")
        else:
            n_ok += 1
            print(f"[same] {rf.name}")
    for rt in sorted(ref.glob("*.txt")):
        new = cohorts_dir / rt.name
        same = new.exists() and new.read_bytes() == rt.read_bytes()
        print(f"[{'same' if same else 'DIFF'}] {rt.name}")
        n_ok += same
        n_bad += not same
    print(f"\n{n_ok} identical, {n_bad} different or missing")
    return 0 if n_bad == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())
