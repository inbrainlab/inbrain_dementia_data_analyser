"""Run context shared by the dataset pipelines."""
from __future__ import annotations

import importlib
import os
import warnings
from dataclasses import dataclass, field
from pathlib import Path
from types import ModuleType
from typing import Callable

import pandas as pd

from ..config import DatasetConfig, ProjectConfig, load_config


@dataclass
class DatasetContext:
    """Everything a dataset pipeline needs to run.

    Attributes
    ----------
    root:
        Root folder that holds one ``Tables_<DATASET>`` folder per dataset.
    dataset:
        Configuration of the dataset being processed.
    config:
        Full project configuration (variables, codebooks, ...).
    files:
        Raw table key -> resolved file (``None`` when the table was not found).
    """

    root: Path
    dataset: DatasetConfig
    config: ProjectConfig
    files: dict[str, Path | None] = field(default_factory=dict)

    @property
    def name(self) -> str:
        return self.dataset.name

    @property
    def folder(self) -> Path:
        """Raw-data folder of the dataset."""
        return self.dataset.folder_path(self.root)

    @property
    def output_dir(self) -> Path:
        return self.dataset.output_dir(self.root)

    @property
    def output_path(self) -> Path:
        return self.dataset.output_path(self.root)

    @property
    def settings(self) -> dict:
        """Cohort-specific settings from the dataset file."""
        return self.dataset.settings

    def path(self, key: str) -> Path | None:
        """Resolved path of a raw table (``None`` if missing)."""
        if key not in self.dataset.tables:
            raise KeyError(f"{self.name}: table {key!r} is not declared in the dataset configuration")
        return self.files.get(key)

    def path_str(self, key: str) -> str:
        """Path of a raw table as text.

        Missing tables return a path that does not exist, so the ported loaders
        report them as skipped exactly as the notebooks did (validation has
        already decided whether that is acceptable).
        """
        p = self.path(key)
        if p is not None:
            return str(p)
        return str(self.folder / f"<missing:{key}>")

    def paths(self) -> dict[str, str]:
        """Every declared raw table -> path text (see :meth:`path_str`)."""
        return {k: self.path_str(k) for k in self.dataset.tables}


def dataset_module(name: str, config: ProjectConfig | None = None) -> ModuleType:
    """Import the pipeline module of a dataset (``datasets/<module>.py``)."""
    cfg = config or load_config()
    ds = cfg.dataset(name)
    return importlib.import_module(f"{__package__}.{ds.module}")


def run_dataset(ctx: DatasetContext) -> pd.DataFrame:
    """Run the pipeline of ``ctx.dataset`` and return the per-subject table."""
    mod = dataset_module(ctx.name, ctx.config)
    run: Callable[[DatasetContext], pd.DataFrame] = getattr(mod, "run")
    os.makedirs(ctx.output_dir, exist_ok=True)
    with warnings.catch_warnings():
        # The ported pipelines rely on pandas behaviours that emit routine
        # warnings (date-format inference, chained assignment, downcasting);
        # the notebooks silenced them globally. Here they are silenced only
        # while a pipeline runs.
        warnings.simplefilter("ignore", category=UserWarning)
        warnings.simplefilter("ignore", category=FutureWarning)
        warnings.simplefilter("ignore", category=pd.errors.SettingWithCopyWarning)
        warnings.simplefilter("ignore", category=pd.errors.PerformanceWarning)
        return run(ctx)
