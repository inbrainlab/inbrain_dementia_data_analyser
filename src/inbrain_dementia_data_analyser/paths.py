"""Locating raw files under the root path.

File names in the source portals carry the download date
(``All_Subjects_DXSUM_08Sep2025.csv``, ``..._2025-03-21.csv``), which changes
at every re-download. Table patterns in the dataset configuration therefore
use two placeholders:

* ``{date}`` - any of the date formats used by the portals
  (``08Sep2025``, ``01-Jun-2018``, ``2025-03-21``, ``20250321``);
* ``*`` - any text.

Matching is case-insensitive and, when the dataset allows it, tolerates the
``Copy of `` prefix that Google Drive adds to copied files. When several files
match, the one with the most recent date in its name wins (file modification
time breaks ties and decides when there is no date).
"""
from __future__ import annotations

import os
import re
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

DATE_RX = (r"(?:\d{1,2}[A-Za-z]{3}\d{4}|\d{1,2}-[A-Za-z]{3}-\d{4}|\d{4}-\d{2}-\d{2}"
           r"|\d{4}\.\d{2}\.\d{2}|\d{8})")
_DATE_FORMATS = ("%d%b%Y", "%d-%b-%Y", "%Y-%m-%d", "%Y.%m.%d", "%Y%m%d")
COPY_PREFIX_RX = r"(?:copy of\s+)?"


def pattern_to_regex(pattern: str, allow_copy_prefix: bool = True) -> re.Pattern:
    """Compile a table pattern (with ``{date}`` and ``*``) to a regular expression."""
    parts = re.split(r"(\{date\}|\*)", pattern)
    rx = "".join(
        f"(?P<date>{DATE_RX})" if p == "{date}" and "(?P<date>" not in "".join(parts[:i])
        else DATE_RX if p == "{date}"
        else ".*" if p == "*"
        else re.escape(p)
        for i, p in enumerate(parts))
    prefix = COPY_PREFIX_RX if allow_copy_prefix else ""
    return re.compile(rf"^{prefix}{rx}$", re.I)


def parse_name_date(text: str | None) -> datetime | None:
    """Parse a portal date token; ``None`` if it is not a date."""
    if not text:
        return None
    for fmt in _DATE_FORMATS:
        try:
            return datetime.strptime(text, fmt)
        except ValueError:
            continue
    return None


@dataclass(frozen=True)
class FileMatch:
    """Result of resolving one table pattern."""

    path: Path | None
    candidates: tuple[Path, ...] = ()

    @property
    def found(self) -> bool:
        return self.path is not None

    @property
    def ambiguous(self) -> bool:
        return len(self.candidates) > 1


def find_files(folder: str | os.PathLike, pattern: str, *, recursive: bool = False,
               allow_copy_prefix: bool = True) -> list[Path]:
    """All files in ``folder`` whose name matches ``pattern``."""
    folder = Path(folder)
    if not folder.is_dir():
        return []
    rx = pattern_to_regex(pattern, allow_copy_prefix)
    it = folder.rglob("*") if recursive else folder.iterdir()
    return sorted(p for p in it if p.is_file() and rx.match(p.name))


def _rank(path: Path, rx: re.Pattern) -> tuple:
    m = rx.match(path.name)
    d = parse_name_date(m.group("date")) if m and "date" in rx.groupindex and m.group("date") else None
    copy = path.name.lower().startswith("copy of ")
    return (d or datetime.min, not copy, path.stat().st_mtime)


def resolve_file(folder: str | os.PathLike, pattern: str, *, recursive: bool = False,
                 allow_copy_prefix: bool = True) -> FileMatch:
    """Pick the file for ``pattern``: newest name date, then non-copy, then mtime."""
    cands = find_files(folder, pattern, recursive=recursive, allow_copy_prefix=allow_copy_prefix)
    if not cands:
        return FileMatch(None, ())
    rx = pattern_to_regex(pattern, allow_copy_prefix)
    best = max(cands, key=lambda p: _rank(p, rx))
    return FileMatch(best, tuple(cands))
