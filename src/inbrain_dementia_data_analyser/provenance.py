"""Provenance records written next to every output."""
from __future__ import annotations

import json
import os
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import pandas as pd

from ._version import __version__
from .config import ProjectConfig
from .io import file_fingerprint


def provenance_path(output: str | os.PathLike) -> Path:
    """``<output>.provenance.json``."""
    p = Path(output)
    return p.with_name(p.name + ".provenance.json")


def write_provenance(output: str | os.PathLike, *, step: str, config: ProjectConfig,
                     inputs: Iterable[str | os.PathLike] = (), table: pd.DataFrame | None = None,
                     parameters: dict[str, Any] | None = None) -> Path:
    """Write a small JSON record describing how ``output`` was produced.

    The record holds the package version, the configuration hash, the input
    files (name, size, modification time), the run parameters and the output
    shape. It contains no subject-level data.
    """
    inputs_info = []
    for p in inputs:
        if p and Path(p).is_file():
            inputs_info.append(file_fingerprint(p))
    rec = {
        "step": step,
        "created": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "package_version": __version__,
        "python": sys.version.split()[0],
        "pandas": pd.__version__,
        "platform": platform.platform(),
        "config_dir": str(config.config_dir),
        "config_sha256": config.fingerprint(),
        "parameters": parameters or {},
        "inputs": inputs_info,
        "output": str(output),
    }
    if table is not None:
        rec["rows"] = int(len(table))
        rec["columns"] = int(table.shape[1])
    path = provenance_path(output)
    path.write_text(json.dumps(rec, indent=2, default=str), encoding="utf-8")
    return path
