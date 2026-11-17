"""Checking the folder structure and raw tables before processing.

A pipeline that cannot find one of its tables still runs and writes an output
that looks complete, with whole domains silently missing. Validation makes
that loud: it resolves every declared table, checks that the identifier and
date columns (and the columns read by configured variables) are in the
header, and returns a report that :func:`~.preprocess.preprocess_dataset`
refuses to ignore unless asked to.
"""
from __future__ import annotations

import logging
import os
from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from .config import ProjectConfig, load_config
from .imaging.ida import ida_find
from .io import TableReadError, read_header
from .paths import resolve_file

log = logging.getLogger(__name__)


@dataclass
class TableCheck:
    """Outcome for one declared table."""

    key: str
    pattern: str
    required: bool
    path: Path | None = None
    candidates: tuple[Path, ...] = ()
    missing_columns: list[str] = field(default_factory=list)
    missing_variable_columns: list[str] = field(default_factory=list)
    error: str | None = None

    @property
    def found(self) -> bool:
        return self.path is not None

    @property
    def status(self) -> str:
        if self.error:
            return "unreadable"
        if not self.found:
            return "missing" if self.required else "optional-missing"
        if self.missing_columns:
            return "missing-columns"
        return "ok"


@dataclass
class ValidationReport:
    """Validation result of one dataset."""

    dataset: str
    root: Path
    folder: Path
    folder_exists: bool
    tables: list[TableCheck] = field(default_factory=list)
    ida_files: list[Path] = field(default_factory=list)
    ida_required: bool = False
    messages: list[str] = field(default_factory=list)

    @property
    def errors(self) -> list[str]:
        """Problems that make the pipeline output incomplete."""
        out = []
        if not self.folder_exists:
            out.append(f"folder not found: {self.folder}")
        for t in self.tables:
            if t.status == "missing":
                out.append(f"required table '{t.key}' not found (pattern {t.pattern!r})")
            elif t.status == "unreadable":
                out.append(f"table '{t.key}' cannot be read: {t.error}")
            elif t.status == "missing-columns":
                out.append(f"table '{t.key}' ({t.path.name}) lacks column(s) {t.missing_columns}")
        if self.ida_required and not self.ida_files:
            out.append("no LONI IDA search export (idaSearch_*.csv with 'Study Date' and "
                       "'Image ID') in the dataset folder")
        return out

    @property
    def warnings(self) -> list[str]:
        """Problems that only affect individual variables or optional tables."""
        out = []
        for t in self.tables:
            if t.status == "optional-missing":
                out.append(f"optional table '{t.key}' not found (pattern {t.pattern!r})")
            if t.missing_variable_columns:
                out.append(f"table '{t.key}': column(s) used by variables.yaml not found "
                           f"{t.missing_variable_columns} (those variables will be skipped)")
            if t.found and len(t.candidates) > 1:
                out.append(f"table '{t.key}': {len(t.candidates)} files match, using "
                           f"{t.path.name}")
        return out + self.messages

    @property
    def ok(self) -> bool:
        return not self.errors

    def files(self) -> dict[str, Path | None]:
        """Table key -> resolved path."""
        return {t.key: t.path for t in self.tables}

    def to_frame(self) -> pd.DataFrame:
        """One row per table (for display or saving)."""
        return pd.DataFrame([{
            "dataset": self.dataset, "table": t.key, "status": t.status,
            "required": t.required, "file": t.path.name if t.path else None,
            "pattern": t.pattern, "n_candidates": len(t.candidates),
            "missing_columns": ", ".join(t.missing_columns),
            "missing_variable_columns": ", ".join(t.missing_variable_columns),
        } for t in self.tables])

    def summary(self) -> str:
        """Human-readable summary."""
        n_ok = sum(t.status == "ok" for t in self.tables)
        lines = [f"{self.dataset}: {'OK' if self.ok else 'FAILED'} - {n_ok}/{len(self.tables)} "
                 f"tables ok; folder {self.folder}"]
        if self.ida_required or self.ida_files:
            lines.append(f"  IDA exports: {[p.name for p in self.ida_files] or 'none'}")
        lines += [f"  ERROR   {e}" for e in self.errors]
        lines += [f"  WARNING {w}" for w in self.warnings]
        return "\n".join(lines)

    def raise_if_failed(self) -> None:
        """Raise :class:`ValidationError` when :attr:`errors` is not empty."""
        if not self.ok:
            raise ValidationError(self)


class ValidationError(RuntimeError):
    """Raised when a dataset fails validation."""

    def __init__(self, report: ValidationReport):
        self.report = report
        super().__init__(report.summary())


def validate_dataset(root_path: str | os.PathLike, dataset: str, *,
                     config: ProjectConfig | None = None, check_columns: bool = True) -> ValidationReport:
    """Validate the raw files of one dataset.

    Parameters
    ----------
    root_path:
        Folder that holds the ``Tables_<DATASET>`` folders.
    dataset:
        Dataset name (``"ADNI"``, ``"PREVENT-AD"``, ...).
    config:
        Project configuration (loaded from the default folder if omitted).
    check_columns:
        Read each table header and check the declared columns.

    Returns
    -------
    ValidationReport
    """
    cfg = config or load_config()
    ds = cfg.dataset(dataset)
    root = Path(root_path).expanduser()
    folder = ds.folder_path(root)
    rep = ValidationReport(dataset=ds.name, root=root, folder=folder, folder_exists=folder.is_dir(),
                           ida_required=bool(ds.ida.get("required", False)))
    if not rep.folder_exists:
        rep.tables = [TableCheck(k, t.pattern, t.required) for k, t in ds.tables.items()]
        return rep
    for key, spec in ds.tables.items():
        m = resolve_file(folder, spec.pattern, recursive=spec.recursive,
                         allow_copy_prefix=ds.allow_copy_prefix)
        chk = TableCheck(key, spec.pattern, spec.required, m.path, m.candidates)
        if m.found and check_columns:
            try:
                header = set(read_header(m.path, sheet=spec.sheet))
            except (TableReadError, OSError, ValueError) as exc:
                chk.error = str(exc)
                header = set()
            if not chk.error:
                lower = {h.lower() for h in header}
                if spec.id_columns and not any(c in header for c in spec.id_columns):
                    chk.missing_columns.append(" or ".join(spec.id_columns))
                if spec.date_columns and not any(c in header for c in spec.date_columns):
                    chk.missing_columns.append(" or ".join(spec.date_columns))
                chk.missing_columns += [c for c in spec.columns if c not in header]
                var_cols = []
                for _, src in cfg.sources(ds.name, key):
                    if src.transform or src.instrument:
                        continue  # derived or wide-timepoint columns are checked when built
                    if not any(c in header or (src.case_insensitive and c.lower() in lower)
                               for c in src.column):
                        var_cols.append(" or ".join(src.column))
                chk.missing_variable_columns = var_cols
        rep.tables.append(chk)
    if ds.uses_ida:
        rep.ida_files = [Path(p) for p in ida_find(str(folder)).get("new", [])]
        if not rep.ida_files and not rep.ida_required:
            rep.messages.append("no LONI IDA search export found; imaging falls back to the "
                                "dataset's own imaging tables")
    return rep


def validate_all(root_path: str | os.PathLike, datasets: list[str] | None = None, *,
                 config: ProjectConfig | None = None) -> dict[str, ValidationReport]:
    """Validate several datasets (all configured datasets by default)."""
    cfg = config or load_config()
    names = datasets or list(cfg.datasets)
    return {cfg.dataset(n).name: validate_dataset(root_path, n, config=cfg) for n in names}


def check_root(root_path: str | os.PathLike, config: ProjectConfig | None = None) -> pd.DataFrame:
    """Which dataset folders exist under ``root_path``."""
    cfg = config or load_config()
    root = Path(root_path).expanduser()
    return pd.DataFrame([{"dataset": d.name, "folder": d.folder,
                          "exists": d.folder_path(root).is_dir()} for d in cfg.datasets.values()])
