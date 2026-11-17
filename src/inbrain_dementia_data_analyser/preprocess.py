"""Stage 1: one integrated per-subject table per dataset."""
from __future__ import annotations

import gc
import logging
import os
import time
from pathlib import Path

import pandas as pd

from .config import ProjectConfig, load_config
from .datasets.base import DatasetContext, run_dataset
from .provenance import write_provenance
from .validation import ValidationError, validate_dataset

log = logging.getLogger(__name__)


def preprocess_dataset(dataset_name: str, root_path: str | os.PathLike, *,
                       config: ProjectConfig | None = None,
                       config_dir: str | os.PathLike | None = None,
                       allow_missing: bool = False) -> pd.DataFrame:
    """Validate and pre-process one dataset.

    Reads the raw tables of ``dataset_name`` under ``root_path`` and writes
    ``<root>/<folder>/<output_folder>/<output_file>`` (for example
    ``Tables_ADNI/Integrated_ADNI/ADNI_longitudinal.csv``) plus a provenance
    record next to it.

    Parameters
    ----------
    dataset_name:
        ``"ADNI"``, ``"A4"``, ``"ABVIB"``, ``"AIBL"``, ``"MCSA"``, ``"NACC"``,
        ``"NIFD"``, ``"OASIS"``, ``"PREVENT-AD"`` or ``"4RTNI"``.
    root_path:
        Folder holding the ``Tables_<DATASET>`` folders.
    config, config_dir:
        Configuration object or folder (default: the repository ``configs/``).
    allow_missing:
        Run even when validation reports missing required tables. The output
        will lack the corresponding domains.

    Returns
    -------
    pandas.DataFrame
        The per-subject table that was written.

    Raises
    ------
    ValidationError
        If validation fails and ``allow_missing`` is False.
    """
    cfg = config or load_config(config_dir)
    ds = cfg.dataset(dataset_name)
    report = validate_dataset(root_path, ds.name, config=cfg)
    for line in report.summary().splitlines():
        (log.warning if ("ERROR" in line or "WARNING" in line) else log.info)(line)
    if not report.ok and not allow_missing:
        raise ValidationError(report)
    ctx = DatasetContext(root=Path(root_path).expanduser(), dataset=ds, config=cfg,
                         files=report.files())
    t0 = time.time()
    master = run_dataset(ctx)
    inputs = [p for p in report.files().values() if p] + list(report.ida_files)
    write_provenance(ctx.output_path, step=f"preprocess:{ds.name}", config=cfg, inputs=inputs,
                     table=master, parameters={"root_path": str(root_path),
                                               "allow_missing": allow_missing,
                                               "seconds": round(time.time() - t0, 1),
                                               "validation_errors": report.errors})
    return master


def preprocess_all(root_path: str | os.PathLike, datasets: list[str] | None = None, *,
                   config: ProjectConfig | None = None,
                   config_dir: str | os.PathLike | None = None,
                   allow_missing: bool = False, stop_on_error: bool = False) -> pd.DataFrame:
    """Pre-process several datasets one after the other.

    The per-dataset tables are written to disk and released from memory after
    each dataset. Returns one summary row per dataset (status, rows, columns,
    output path, error message).
    """
    cfg = config or load_config(config_dir)
    names = datasets or list(cfg.datasets)
    rows = []
    for name in names:
        ds = cfg.dataset(name)
        try:
            df = preprocess_dataset(ds.name, root_path, config=cfg, allow_missing=allow_missing)
            rows.append({"dataset": ds.name, "status": "ok", "rows": len(df), "columns": df.shape[1],
                         "output": str(ds.output_path(root_path)), "error": ""})
            del df
        except (ValidationError, FileNotFoundError, KeyError, ValueError, OSError) as exc:
            log.error("[ERROR] %s failed: %s", ds.name, exc)
            rows.append({"dataset": ds.name, "status": "failed", "rows": 0, "columns": 0,
                         "output": "", "error": str(exc).splitlines()[0]})
            if stop_on_error:
                raise
        gc.collect()
    return pd.DataFrame(rows)
