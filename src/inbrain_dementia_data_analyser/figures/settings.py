"""Figure settings: ``configs/figures.yaml`` merged with a user spec."""
from __future__ import annotations

import copy
import os
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from ..config import ConfigError, ProjectConfig, load_config

#: used when the configuration folder has no figures.yaml (older copies)
BUILTIN_DEFAULTS: dict[str, Any] = {
    "diagnosis_groups": {
        "order": ["CN", "MCI", "AD", "LBD", "FTD / tauopathy", "VaD", "Mixed"],
        "colors": {"CN": "#2a78d6", "MCI": "#eb6834", "AD": "#1baf7a", "LBD": "#eda100",
                   "FTD / tauopathy": "#e87ba4", "VaD": "#008300", "Mixed": "#6250d6"},
        "labels": {"CN": "CN", "CON": "CN", "MCI": "MCI", "AD": "AD", "VaD": "VaD", "LBD": "LBD",
                   "LBD+PD": "LBD", "bvFTD": "FTD / tauopathy", "FTD": "FTD / tauopathy",
                   "svPPA": "FTD / tauopathy", "nfvPPA": "FTD / tauopathy", "PSP": "FTD / tauopathy",
                   "CBS": "FTD / tauopathy", "FTD+PSP": "FTD / tauopathy"},
        "compound_group": "Mixed",
        "stage_rank": {"CN": 0, "MCI": 1},
        "default_rank": 2,
        "light_text": ["CN", "MCI", "AD", "VaD", "Mixed"],
    },
    "modalities": {"T1w": {"measure": "MRI_T1w"}, "CDR": {"measure": "CDR"},
                   "MMSE": {"measure": "MMSE"}, "MoCA": {"measure": "MoCA"},
                   "p-tau": {"measure": "CSF_PTAU"}, "Aβ42": {"measure": "CSF_AB42"}},
    "reference": "T1w",
    "bundles": {"p-tau + Aβ42": ["p-tau", "Aβ42"]},
    "composites": "auto",
    "venn": ["MMSE", "MoCA"],
    "longitudinal": {"window_months": 6, "max_months": 60,
                     "availability": ["T1w", "CDR", "MMSE", "MoCA", "p-tau"],
                     "sankey_min_flow": 3, "upset_min_n": 5},
    "style": {"dpi": 300, "formats": ["png"], "font": "DejaVu Sans"},
}

#: keys whose value replaces the default as a whole (no deep merge)
_REPLACE = {"modalities", "bundles", "composites", "venn", "labels", "order", "colors",
            "availability", "formats", "light_text", "stage_rank"}


def _merge(base: dict, over: dict) -> dict:
    out = copy.deepcopy(base)
    for k, v in (over or {}).items():
        if k in out and isinstance(out[k], dict) and isinstance(v, dict) and k not in _REPLACE:
            out[k] = _merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


def load_figure_settings(spec: str | os.PathLike | dict | None = None,
                         config: ProjectConfig | None = None) -> dict[str, Any]:
    """Defaults of ``figures.yaml`` (configuration folder) updated with ``spec``.

    ``spec`` is a YAML file or a dict with any of the keys of ``figures.yaml``;
    ``modalities``, ``bundles``, ``venn``, ``composites``, colour and label maps
    replace the defaults as a whole, the other sections are merged key by key.
    """
    cfg = config or load_config()
    path = Path(cfg.config_dir) / "figures.yaml"
    base = BUILTIN_DEFAULTS
    if path.is_file():
        base = _merge(BUILTIN_DEFAULTS, yaml.safe_load(path.read_text(encoding="utf-8")) or {})
    if spec is None:
        over = {}
    elif isinstance(spec, dict):
        over = spec
    else:
        over = yaml.safe_load(Path(spec).read_text(encoding="utf-8")) or {}
    out = _merge(base, over)
    out["integration_window"] = int(cfg.integration.get("expand", {}).get("window_step_months", 6))
    _validate(out)
    return out


def _validate(s: dict) -> None:
    g = s["diagnosis_groups"]
    missing = [x for x in g["order"] if x not in g["colors"]]
    if missing:
        raise ConfigError(f"figures: diagnosis_groups.colors has no colour for {missing}")
    bad = sorted({v for v in g["labels"].values() if v not in g["order"]})
    if bad:
        raise ConfigError(f"figures: diagnosis_groups.labels maps to groups not in 'order': {bad}")
    mods = s["modalities"] or {}
    if not mods:
        raise ConfigError("figures: 'modalities' is empty")
    for lab, m in mods.items():
        if isinstance(m, str):
            mods[lab] = {"measure": m}
        elif not isinstance(m, dict):
            raise ConfigError(f"figures: modalities.{lab} must be a mapping (measure / column / file)")
    ref = s.get("reference")
    if ref is not None and ref not in mods:
        raise ConfigError(f"figures: reference {ref!r} is not one of the modalities {list(mods)}")
    for b, members in (s.get("bundles") or {}).items():
        unknown = [m for m in members if m not in mods]
        if unknown:
            raise ConfigError(f"figures: bundle {b!r} uses unknown modalities {unknown}")
    v = s.get("venn")
    if v is not None and (len(v) != 2 or any(x not in mods for x in v)):
        raise ConfigError(f"figures: venn must name two modalities of {list(mods)}, got {v}")
    for lab in (s["longitudinal"].get("availability") or []):
        if lab not in mods:
            raise ConfigError(f"figures: longitudinal.availability uses unknown modality {lab!r}")
    w = int(s["longitudinal"]["window_months"])
    if w <= 0:
        raise ConfigError("figures: longitudinal.window_months must be > 0")


class DxGrouper:
    """Canonical diagnosis label -> diagnostic group of the figures."""

    def __init__(self, settings: dict):
        g = settings["diagnosis_groups"]
        self.order: list[str] = list(g["order"])
        self.colors: dict[str, str] = dict(g["colors"])
        self.labels: dict[str, str] = {str(k): v for k, v in g["labels"].items()}
        self._upper = {k.upper(): v for k, v in self.labels.items()}
        self.compound = g.get("compound_group")
        self.rank = dict(g.get("stage_rank") or {})
        self.default_rank = g.get("default_rank", 2)
        self.light_text = set(g.get("light_text") or [])

    def __call__(self, label: Any) -> str | None:
        if label is None or (isinstance(label, float) and pd.isna(label)):
            return None
        s = str(label).strip()
        if not s or s.lower() in ("nan", "none", "nat", "na"):
            return None
        if s in self.labels:
            return self.labels[s]
        if s.upper() in self._upper:
            return self._upper[s.upper()]
        if "+" in s and self.compound:
            return self.compound
        return None

    def rank_of(self, group: str) -> int:
        return self.rank.get(group, self.default_rank)

    def transition(self, a: str, b: str) -> str:
        d = self.rank_of(b) - self.rank_of(a)
        return "Progression" if d > 0 else "Reversion" if d < 0 else "Change of aetiology"
