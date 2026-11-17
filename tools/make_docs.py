#!/usr/bin/env python
"""Write docs/variables.md and docs/datasets.md from configs/ (run after editing the YAML files)."""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO / "src"))
from inbrain_dementia_data_analyser.config import load_config  # noqa: E402


def source_text(cfg, ds, s) -> str:
    table = s.table or cfg.datasets[ds].default_table
    if s.instrument:
        what = f"`{s.instrument}_<tp>_{s.suffix}`" if s.suffix else f"`{s.instrument}_<tp>`"
        if s.transform:
            what = f"{s.transform} of `{s.instrument}_<tp>_<item>` ({len(s.items)} items)"
    elif s.transform:
        n = len(s.items)
        what = f"{s.transform}({', '.join(f'`{i}`' for i in s.items[:4])}{', ...' if n > 4 else ''})"
        if s.transform.startswith("bmi"):
            what = f"{s.transform} of `{s.params.get('weight')}` / `{s.params.get('height')}`"
    else:
        what = " or ".join(f"`{c}`" for c in s.column)
    extra = []
    if s.date_column:
        extra.append(f"date `{s.date_column}`")
    if s.valid_codes is not None:
        extra.append(f"codes {list(s.valid_codes)}")
    if s.missing_values:
        extra.append("sentinels removed")
    if s.row_filter:
        extra.append(f"rows where `{s.row_filter['column']}` ~ `{s.row_filter['pattern']}`")
    return f"{table}: {what}" + (f" ({'; '.join(extra)})" if extra else "")


def main() -> None:
    cfg = load_config(REPO / "configs")
    order = list(cfg.datasets)
    lines = ["# Variables", "",
             "Generated from `configs/variables.yaml` by `tools/make_variable_docs.py`; do not edit by hand.",
             "", f"{len(cfg.variables)} variables. *expand* = placed into time windows in the expanded master.",
             "", "| Variable | Label | Domain | Expand | Sources |", "|---|---|---|---|---|"]
    for name, v in cfg.variables.items():
        srcs = "<br>".join(f"**{ds}** {source_text(cfg, ds, s)}" for ds in order
                           for s in v.sources.get(ds, ()))
        if v.builtin:
            srcs = (srcs + "<br>" if srcs else "") + "*built by dataset code*"
        lines.append(f"| `{name}` | {v.label} | {v.domain} | {'yes' if v.expand else ''} | {srcs} |")
    lines += ["", "## Imaging types", "", "| Type | Modality | Description | Expand |", "|---|---|---|---|"]
    for t in cfg.imaging["types"]:
        lines.append(f"| `{t['name']}` | {t.get('modality', '')} | {t.get('description', '')} | "
                     f"{'yes' if t.get('expand', True) else ''} |")
    (REPO / "docs" / "variables.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"wrote {REPO / 'docs' / 'variables.md'}")

    d = ["# Datasets and raw files", "",
         "Generated from `configs/datasets/*.yaml` by `tools/make_docs.py`; do not edit by hand.", "",
         "File patterns: `{date}` = any download date (`08Sep2025`, `01-Jun-2018`, `2025-03-21`, "
         "`20250321`); `*` = any text. Matching ignores case and accepts the `Copy of ` prefix "
         "added by Google Drive. When several files match, the newest date in the name wins.", ""]
    for ds in cfg.datasets.values():
        d += [f"## {ds.name} - {ds.label}", "",
              f"* Portal: {ds.portal}" if ds.portal else "",
              f"* Raw folder: `<root>/{ds.folder}/`",
              f"* Output: `<root>/{ds.folder}/{ds.output_folder}/{ds.output_file}`",
              f"* Subject ID: `{ds.subject_id}`; time axis: {ds.time_axis}",
              (f"* LONI IDA search export: `idaSearch_*.csv` in the raw folder "
               f"({'required' if ds.ida.get('required') else 'optional; the dataset tables are the fallback'})"
               if ds.uses_ida else ""), "",
              "| Table | File pattern | Required | ID column | Date column | Notes |",
              "|---|---|---|---|---|---|"]
        for k, t in ds.tables.items():
            d.append(f"| `{k}` | `{t.pattern}` | {'yes' if t.required else 'no'} | "
                     f"{' / '.join(t.id_columns)} | {' / '.join(t.date_columns)} | {t.description} |")
        d.append("")
    (REPO / "docs" / "datasets.md").write_text("\n".join(x for x in d if x is not None) + "\n", encoding="utf-8")
    print(f"wrote {REPO / 'docs' / 'datasets.md'}")


if __name__ == "__main__":
    main()
