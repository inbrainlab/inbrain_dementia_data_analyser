"""Configuration loading and validation.

All harmonisation choices live in YAML files under the ``configs/`` folder of
the repository (or any folder given by ``config_dir`` / the environment
variable ``INBRAIN_CONFIG_DIR``):

* ``datasets/<name>.yaml`` - folder names, raw tables and the columns that must
  be present, subject ID, time axis and cohort-specific settings;
* ``variables.yaml`` - every scored variable (cognitive test, scale, vital
  sign, lab value) with the column that holds it in each dataset;
* ``harmonization.yaml`` - column renames applied when the cohorts are merged;
* ``demographics.yaml`` - sex / race / ethnicity / handedness / APOE codebooks;
* ``diagnosis.yaml`` - raw diagnosis label -> canonical label;
* ``imaging.yaml`` - canonical imaging types and measure groups;
* ``integration.yaml`` - settings of the time-window expansion.

The files are validated when they are loaded; a problem raises
:class:`ConfigError` with the file, the entry and what is wrong.
"""
from __future__ import annotations

import hashlib
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

import yaml

ENV_CONFIG_DIR = "INBRAIN_CONFIG_DIR"

SOURCE_KEYS = {
    "table", "column", "date_column", "transform", "items", "params",
    "missing_values", "valid_range", "range_inclusive", "valid_codes",
    "row_filter", "case_insensitive", "instrument", "suffix", "date_suffix",
    "note",
}
VARIABLE_KEYS = {
    "label", "domain", "kind", "unit", "expand", "builtin", "description",
    "sources",
}
DATASET_KEYS = {
    "name", "label", "portal", "description", "folder", "output_folder",
    "output_file", "subject_id", "time_axis", "module", "default_table",
    "allow_copy_prefix", "ida", "tables", "settings",
}
TABLE_KEYS = {
    "pattern", "required", "recursive", "id_columns", "date_columns",
    "columns", "sheet", "description", "generic", "date_column", "id_column",
}
TIME_AXES = {"date", "days", "months", "visit"}
VARIABLE_KINDS = {"score", "categorical", "flag"}


class ConfigError(ValueError):
    """Raised when a configuration file is missing or invalid."""


# ---------------------------------------------------------------------------
# data classes
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class TableSpec:
    """One raw input table of a dataset."""

    key: str
    pattern: str
    required: bool = True
    recursive: bool = False
    id_columns: tuple[str, ...] = ()
    date_columns: tuple[str, ...] = ()
    columns: tuple[str, ...] = ()
    sheet: str | int | None = None
    description: str = ""
    generic: bool = False
    id_column: str | None = None
    date_column: str | None = None


@dataclass(frozen=True)
class DatasetConfig:
    """Folder layout, raw tables and settings of one dataset."""

    name: str
    folder: str
    output_folder: str
    output_file: str
    subject_id: str
    time_axis: str
    module: str
    label: str = ""
    portal: str = ""
    description: str = ""
    default_table: str | None = None
    allow_copy_prefix: bool = True
    ida: dict[str, Any] = field(default_factory=dict)
    tables: dict[str, TableSpec] = field(default_factory=dict)
    settings: dict[str, Any] = field(default_factory=dict)

    @property
    def uses_ida(self) -> bool:
        """True when the dataset reads a LONI IDA search export."""
        return bool(self.ida.get("enabled", False))

    def folder_path(self, root: str | os.PathLike) -> Path:
        """Raw-data folder of the dataset under ``root``."""
        return Path(root) / self.folder

    def output_dir(self, root: str | os.PathLike) -> Path:
        """Folder that receives the per-dataset integrated table."""
        return self.folder_path(root) / self.output_folder

    def output_path(self, root: str | os.PathLike) -> Path:
        """Full path of the per-dataset integrated table."""
        return self.output_dir(root) / self.output_file


@dataclass(frozen=True)
class SourceSpec:
    """Where one variable lives in one dataset, and how to clean it."""

    variable: str
    dataset: str
    table: str | None
    column: tuple[str, ...] = ()
    date_column: str | None = None
    transform: str | None = None
    items: tuple[str, ...] = ()
    params: dict[str, Any] = field(default_factory=dict)
    missing_values: tuple[Any, ...] = ()
    valid_range: tuple[float | None, float | None] | None = None
    range_inclusive: bool = True
    valid_codes: tuple[Any, ...] | None = None
    row_filter: dict[str, Any] | None = None
    case_insensitive: bool = False
    instrument: str | None = None
    suffix: str | None = None
    date_suffix: str | None = None

    @property
    def needs_preparation(self) -> bool:
        """True when the raw column cannot be used as it is."""
        return bool(self.transform or self.missing_values or self.valid_range
                    or self.valid_codes is not None or self.row_filter)

    def raw_columns(self) -> list[str]:
        """Every raw column this source may read (for header checks)."""
        cols = list(self.column) + list(self.items)
        for key in ("weight", "height", "weight_unit_column",
                    "height_unit_column", "bmi_column", "numerator",
                    "denominator", "filter_column"):
            val = self.params.get(key)
            if isinstance(val, str):
                cols.append(val)
            elif isinstance(val, (list, tuple)):
                cols.extend(str(v) for v in val)
        if self.row_filter and self.row_filter.get("column"):
            cols.append(self.row_filter["column"])
        return cols


@dataclass(frozen=True)
class VariableConfig:
    """One harmonised output variable."""

    name: str
    label: str = ""
    domain: str = ""
    kind: str = "score"
    unit: str = ""
    expand: bool = False
    builtin: bool = False
    description: str = ""
    sources: dict[str, tuple[SourceSpec, ...]] = field(default_factory=dict)


@dataclass
class ProjectConfig:
    """Everything loaded from a configuration folder."""

    config_dir: Path
    datasets: dict[str, DatasetConfig]
    variables: dict[str, VariableConfig]
    harmonization: dict[str, Any]
    demographics: dict[str, Any]
    diagnosis: dict[str, Any]
    imaging: dict[str, Any]
    integration: dict[str, Any]
    files: list[Path] = field(default_factory=list)

    # -- lookups -----------------------------------------------------------
    def dataset(self, name: str) -> DatasetConfig:
        """Return a dataset configuration (case-insensitive name, aliases allowed)."""
        key = canonical_dataset_name(name, self.datasets)
        return self.datasets[key]

    def sources(self, dataset: str, table: str | None = None,
                *, include_builtin: bool = False) -> list[tuple[str, SourceSpec]]:
        """All ``(variable_name, source)`` pairs of ``dataset``.

        Parameters
        ----------
        dataset:
            Dataset name.
        table:
            Restrict to one table key. ``None`` returns every source.
        include_builtin:
            Also return sources of variables flagged ``builtin`` (those are
            produced by dataset code and listed only for documentation).

        Returns
        -------
        list of tuple
            In the order the variables appear in ``variables.yaml``.
        """
        ds = canonical_dataset_name(dataset, self.datasets)
        out = []
        for name, var in self.variables.items():
            if var.builtin and not include_builtin:
                continue
            for src in var.sources.get(ds, ()):
                eff_table = src.table or self.datasets[ds].default_table
                if table is None or eff_table == table:
                    out.append((name, src))
        return out

    def source_columns(self, dataset: str, table: str | None = None) -> list[str]:
        """Raw columns read by the configured variables of a dataset/table."""
        cols: list[str] = []
        for _, src in self.sources(dataset, table):
            cols.extend(src.raw_columns())
        return list(dict.fromkeys(cols))

    def expand_features(self) -> dict[str, str]:
        """Measures expanded into time windows: ``{stem: 'scored'|'imaging'}``."""
        feats = {name: "scored" for name, v in self.variables.items() if v.expand}
        for t in self.imaging.get("types", []):
            if t.get("expand", True):
                feats[t["name"]] = "imaging"
        return feats

    def fingerprint(self) -> str:
        """SHA-256 of every configuration file (for provenance records)."""
        h = hashlib.sha256()
        for p in sorted(self.files):
            h.update(p.name.encode())
            h.update(p.read_bytes())
        return h.hexdigest()


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------
DATASET_ALIASES = {
    "PREVENT-AD": "PREVENTAD", "PREVENT_AD": "PREVENTAD", "4TRNI": "4RTNI",
    "OASIS3": "OASIS", "OASIS-3": "OASIS", "RTNI4": "4RTNI",
}


def canonical_dataset_name(name: str, known: Iterable[str]) -> str:
    """Map a user-supplied dataset name to its configured key.

    Raises
    ------
    ConfigError
        If the name does not match any configured dataset.
    """
    known = list(known)
    up = str(name).strip().upper()
    up = DATASET_ALIASES.get(up, up)
    for k in known:
        if k.upper() == up:
            return k
    raise ConfigError(f"Unknown dataset {name!r}. Configured datasets: {', '.join(known)}")


def default_config_dir() -> Path:
    """Locate the configuration folder.

    Order: environment variable ``INBRAIN_CONFIG_DIR``; a ``configs`` folder
    next to the source tree (clone / editable install); the copy installed
    under ``<prefix>/share/inbrain_dementia_data_analyser/configs``; a
    ``configs`` folder in the working directory.
    """
    env = os.environ.get(ENV_CONFIG_DIR)
    if env:
        return Path(env).expanduser()
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "configs"
        if (cand / "variables.yaml").is_file():
            return cand
    import sys

    for base in (Path(sys.prefix), Path(sys.base_prefix), Path.home() / ".local"):
        cand = base / "share" / "inbrain_dementia_data_analyser" / "configs"
        if (cand / "variables.yaml").is_file():
            return cand
    cand = Path.cwd() / "configs"
    if (cand / "variables.yaml").is_file():
        return cand
    raise ConfigError(
        "Could not find the configuration folder. Install the package from a "
        "clone with `pip install -e .`, pass config_dir=..., or set "
        f"{ENV_CONFIG_DIR}.")


def _read_yaml(path: Path) -> Any:
    if not path.is_file():
        raise ConfigError(f"Configuration file not found: {path}")
    try:
        with path.open(encoding="utf-8") as fh:
            return yaml.safe_load(fh) or {}
    except yaml.YAMLError as exc:
        raise ConfigError(f"{path.name}: invalid YAML ({exc})") from exc


def _as_tuple(val: Any) -> tuple:
    if val is None:
        return ()
    if isinstance(val, (list, tuple)):
        return tuple(val)
    return (val,)


def _check_keys(where: str, data: dict, allowed: set[str]) -> None:
    unknown = set(data) - allowed
    if unknown:
        raise ConfigError(f"{where}: unknown key(s) {sorted(unknown)}; allowed: {sorted(allowed)}")


def _parse_table(ds: str, key: str, raw: Any, fname: str) -> TableSpec:
    where = f"{fname} -> tables.{key}"
    if isinstance(raw, str):
        raw = {"pattern": raw}
    if not isinstance(raw, dict):
        raise ConfigError(f"{where}: expected a mapping or a file pattern, got {type(raw).__name__}")
    _check_keys(where, raw, TABLE_KEYS)
    if not raw.get("pattern"):
        raise ConfigError(f"{where}: 'pattern' is required (file name with optional {{date}} and *)")
    generic = bool(raw.get("generic", False))
    id_col = raw.get("id_column")
    date_col = raw.get("date_column")
    if generic and not (id_col and date_col):
        raise ConfigError(f"{where}: generic tables need 'id_column' and 'date_column'")
    id_cols = _as_tuple(raw.get("id_columns")) or ((id_col,) if id_col else ())
    date_cols = _as_tuple(raw.get("date_columns")) or ((date_col,) if date_col else ())
    return TableSpec(
        key=key, pattern=str(raw["pattern"]), required=bool(raw.get("required", True)),
        recursive=bool(raw.get("recursive", False)),
        id_columns=tuple(map(str, id_cols)), date_columns=tuple(map(str, date_cols)),
        columns=tuple(map(str, _as_tuple(raw.get("columns")))), sheet=raw.get("sheet"),
        description=str(raw.get("description", "")), generic=generic,
        id_column=id_col, date_column=date_col)


def _parse_dataset(path: Path) -> DatasetConfig:
    raw = _read_yaml(path)
    fname = f"datasets/{path.name}"
    if not isinstance(raw, dict):
        raise ConfigError(f"{fname}: expected a mapping at top level")
    _check_keys(fname, raw, DATASET_KEYS)
    for req in ("name", "folder", "output_folder", "output_file", "subject_id", "time_axis", "module"):
        if not raw.get(req):
            raise ConfigError(f"{fname}: '{req}' is required")
    if raw["time_axis"] not in TIME_AXES:
        raise ConfigError(f"{fname}: time_axis must be one of {sorted(TIME_AXES)}, got {raw['time_axis']!r}")
    tables = {str(k): _parse_table(raw["name"], str(k), v, fname)
              for k, v in (raw.get("tables") or {}).items()}
    default_table = raw.get("default_table")
    if default_table and default_table not in tables:
        raise ConfigError(f"{fname}: default_table {default_table!r} is not a declared table")
    return DatasetConfig(
        name=str(raw["name"]), folder=str(raw["folder"]), output_folder=str(raw["output_folder"]),
        output_file=str(raw["output_file"]), subject_id=str(raw["subject_id"]),
        time_axis=str(raw["time_axis"]), module=str(raw["module"]), label=str(raw.get("label", "")),
        portal=str(raw.get("portal", "")), description=str(raw.get("description", "")),
        default_table=default_table, allow_copy_prefix=bool(raw.get("allow_copy_prefix", True)),
        ida=dict(raw.get("ida") or {}), tables=tables, settings=dict(raw.get("settings") or {}))


def _parse_source(var: str, ds: str, raw: Any, datasets: dict[str, DatasetConfig],
                  transforms: set[str]) -> SourceSpec:
    where = f"variables.yaml -> {var}.sources.{ds}"
    if isinstance(raw, str):
        raw = {"column": raw}
    if not isinstance(raw, dict):
        raise ConfigError(f"{where}: expected a mapping, got {type(raw).__name__}")
    _check_keys(where, raw, SOURCE_KEYS)
    dcfg = datasets[ds]
    table = raw.get("table")
    eff = table or dcfg.default_table
    if eff is None:
        raise ConfigError(f"{where}: 'table' is required ({ds} has no default_table). "
                          f"Declared tables: {sorted(dcfg.tables)}")
    if eff not in dcfg.tables:
        raise ConfigError(f"{where}: table {eff!r} is not declared in datasets/{ds.lower()}.yaml "
                          f"(declared: {sorted(dcfg.tables)})")
    transform = raw.get("transform")
    if transform is not None and transform not in transforms:
        raise ConfigError(f"{where}: unknown transform {transform!r}; registered: {sorted(transforms)}")
    column = _as_tuple(raw.get("column"))
    items = _as_tuple(raw.get("items"))
    if not column and not transform and not raw.get("instrument"):
        raise ConfigError(f"{where}: give 'column' (or 'instrument' for wide timepoint tables) "
                          "or a 'transform'")
    vr = raw.get("valid_range")
    if vr is not None:
        if not (isinstance(vr, (list, tuple)) and len(vr) == 2):
            raise ConfigError(f"{where}: valid_range must be [min, max] (use null for open ends)")
        vr = (vr[0], vr[1])
    rf = raw.get("row_filter")
    if rf is not None and not (isinstance(rf, dict) and rf.get("column") and rf.get("pattern")):
        raise ConfigError(f"{where}: row_filter needs 'column' and 'pattern' (regular expression)")
    vc = raw.get("valid_codes")
    return SourceSpec(
        variable=var, dataset=ds, table=table, column=tuple(map(str, column)),
        date_column=raw.get("date_column"), transform=transform,
        items=tuple(map(str, items)), params=dict(raw.get("params") or {}),
        missing_values=_as_tuple(raw.get("missing_values")), valid_range=vr,
        range_inclusive=bool(raw.get("range_inclusive", True)),
        valid_codes=None if vc is None else _as_tuple(vc), row_filter=rf,
        case_insensitive=bool(raw.get("case_insensitive", False)),
        instrument=raw.get("instrument"),
        suffix=None if raw.get("suffix") is None else str(raw.get("suffix")),
        date_suffix=raw.get("date_suffix"))


def _parse_variables(path: Path, datasets: dict[str, DatasetConfig]) -> dict[str, VariableConfig]:
    from .transforms import TRANSFORMS  # local import: transforms imports nothing from here

    raw = _read_yaml(path)
    if not isinstance(raw, dict) or "variables" not in raw:
        raise ConfigError("variables.yaml: expected a top-level 'variables:' mapping")
    out: dict[str, VariableConfig] = {}
    for name, body in (raw["variables"] or {}).items():
        name = str(name)
        where = f"variables.yaml -> {name}"
        if not name.replace("_", "").isalnum():
            raise ConfigError(f"{where}: variable names may only use letters, digits and '_'")
        body = body or {}
        if not isinstance(body, dict):
            raise ConfigError(f"{where}: expected a mapping")
        _check_keys(where, body, VARIABLE_KEYS)
        kind = body.get("kind", "score")
        if kind not in VARIABLE_KINDS:
            raise ConfigError(f"{where}: kind must be one of {sorted(VARIABLE_KINDS)}")
        sources: dict[str, tuple[SourceSpec, ...]] = {}
        for ds, sraw in (body.get("sources") or {}).items():
            try:
                ds_key = canonical_dataset_name(ds, datasets)
            except ConfigError as exc:
                raise ConfigError(f"{where}.sources: {exc}") from None
            items = sraw if isinstance(sraw, list) else [sraw]
            sources[ds_key] = tuple(_parse_source(name, ds_key, s, datasets, set(TRANSFORMS))
                                    for s in items)
        if not sources and not body.get("builtin"):
            raise ConfigError(f"{where}: no sources given (add at least one dataset, "
                              "or set builtin: true for measures produced by dataset code)")
        out[name] = VariableConfig(
            name=name, label=str(body.get("label", "")), domain=str(body.get("domain", "")),
            kind=kind, unit=str(body.get("unit", "")), expand=bool(body.get("expand", False)),
            builtin=bool(body.get("builtin", False)),
            description=str(body.get("description", "")), sources=sources)
    return out


_CACHE: dict[Path, ProjectConfig] = {}


def load_config(config_dir: str | os.PathLike | None = None, *, refresh: bool = False) -> ProjectConfig:
    """Load and validate every configuration file.

    Parameters
    ----------
    config_dir:
        Folder holding the YAML files. Defaults to :func:`default_config_dir`.
    refresh:
        Re-read the files even if this folder was loaded before.

    Returns
    -------
    ProjectConfig

    Raises
    ------
    ConfigError
        If a file is missing or an entry is invalid.
    """
    cdir = Path(config_dir).expanduser().resolve() if config_dir else default_config_dir().resolve()
    if not refresh and cdir in _CACHE:
        return _CACHE[cdir]
    if not cdir.is_dir():
        raise ConfigError(f"Configuration folder not found: {cdir}")
    files: list[Path] = []
    ds_dir = cdir / "datasets"
    datasets: dict[str, DatasetConfig] = {}
    for p in sorted(ds_dir.glob("*.yaml")):
        d = _parse_dataset(p)
        if d.name in datasets:
            raise ConfigError(f"datasets/{p.name}: dataset name {d.name!r} is defined twice")
        datasets[d.name] = d
        files.append(p)
    if not datasets:
        raise ConfigError(f"No dataset files found in {ds_dir}")
    variables = _parse_variables(cdir / "variables.yaml", datasets)
    files.append(cdir / "variables.yaml")
    other = {}
    for name in ("harmonization", "demographics", "diagnosis", "imaging", "integration"):
        p = cdir / f"{name}.yaml"
        other[name] = _read_yaml(p)
        files.append(p)
    cfg = ProjectConfig(config_dir=cdir, datasets=datasets, variables=variables, files=files, **other)
    _validate_cross(cfg)
    _CACHE[cdir] = cfg
    return cfg


def _validate_cross(cfg: ProjectConfig) -> None:
    """Checks that involve several files."""
    types = cfg.imaging.get("types")
    if not isinstance(types, list) or not all(isinstance(t, dict) and t.get("name") for t in types):
        raise ConfigError("imaging.yaml: 'types' must be a list of {name: ...} entries")
    names = [t["name"] for t in types]
    clash = set(names) & {n for n, v in cfg.variables.items() if v.expand}
    if clash:
        raise ConfigError(f"imaging.yaml / variables.yaml: names used both as imaging type and "
                          f"scored variable: {sorted(clash)}")
    for key in ("family_rename", "column_rename"):
        block = cfg.harmonization.get(key, {}) or {}
        for ds in block:
            canonical_dataset_name(ds, cfg.datasets)
