"""Result container of the cohort tools."""
from __future__ import annotations

from dataclasses import dataclass, field

import pandas as pd


@dataclass
class CohortResult:
    """What a cohort tool produced.

    Attributes
    ----------
    table:
        The cohort table (one row per subject, event or visit).
    images:
        Long list of the matched images (one row per subject x measure x image).
    summary:
        Compact per-event sheet (anchor pairing only).
    files:
        Output name -> path of every file written.
    """

    table: pd.DataFrame
    images: pd.DataFrame = field(default_factory=pd.DataFrame)
    summary: pd.DataFrame | None = None
    files: dict[str, str] = field(default_factory=dict)

    @property
    def n_subjects(self) -> int:
        """Number of distinct subjects in :attr:`table`."""
        if "uid" in self.table.columns:
            return int(self.table["uid"].nunique())
        return int(len(self.table))

    def image_ids(self) -> list[str]:
        """Unique image IDs (to paste into the LONI IDA search)."""
        col = next((c for c in ("loniID", "image_id") if c in self.images.columns), None)
        if col is None:
            return []
        return self.images[col].astype(str).str.strip().drop_duplicates().tolist()
