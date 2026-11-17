"""Stage 2: merge the per-dataset tables and expand them into time windows."""
from __future__ import annotations

import gc
import logging
import os
from pathlib import Path

from ..config import ProjectConfig, load_config
from .expand import expand_master
from .merge import merge_cohorts

log = logging.getLogger(__name__)

__all__ = ["integrate", "merge_cohorts", "expand_master", "integrated_paths"]


def integrated_paths(root_path: str | os.PathLike, config: ProjectConfig | None = None) -> dict[str, Path]:
    """Paths of the merged and expanded master tables under ``root_path``."""
    cfg = config or load_config()
    out = Path(root_path) / cfg.integration.get("output_folder", "Integrated_ALL")
    return {"folder": out, "master": out / "MASTER_all_cohorts.csv",
            "coverage": out / "MASTER_column_coverage.csv",
            "expanded": out / "MASTER_expanded_240m.csv"}


def integrate(root_path: str | os.PathLike, datasets: list[str] | None = None, *,
              config: ProjectConfig | None = None, config_dir: str | os.PathLike | None = None,
              expand: bool = True) -> dict[str, Path]:
    """Merge the per-dataset tables, then build the expanded master.

    Parameters
    ----------
    root_path:
        Root folder of the datasets.
    datasets:
        Datasets to merge (default: every configured dataset whose
        pre-processed table exists).
    config, config_dir:
        Configuration object or folder.
    expand:
        Also build ``MASTER_expanded_240m.csv`` (needed by the cohort tools).

    Returns
    -------
    dict
        ``master``, ``coverage``, ``expanded`` and ``folder`` paths.
    """
    cfg = config or load_config(config_dir)
    paths = integrated_paths(root_path, cfg)
    master = merge_cohorts(root_path, datasets, out_dir=paths["folder"], config=cfg)
    log.info("merged: %d subjects x %d columns", len(master), master.shape[1])
    del master
    gc.collect()
    if expand:
        expand_master(paths["master"], paths["expanded"], config=cfg)
    return paths
