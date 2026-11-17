"""Session-aligned image IDs for builders that read per-image tables
(ADNI, OASIS, NACC, PREVENT-AD): one representative ID per distinct date,
plus every ID of the session.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import logging as _logging
from ..logging_utils import make_print as _make_print

log = _logging.getLogger(__name__)
_print = _make_print(log)

# =============================================================================
# SESSION-ALIGNED IMAGE IDs  (ADNI / OASIS / NACC / PREVENT-AD builders)
# =============================================================================
def session_align_ids(dates, ids):
    """
    dates, ids: one entry per IMAGE ROW of one subject x type (same length,
    any order). Returns three lists with ONE entry per distinct date, in the
    order the dates first appear:
        session dates · representative Image ID · all Image IDs ('/'-joined)
    The representative is the lowest numeric ID of the date (first ID when the
    IDs are not numeric); 'NA' when the date has no ID.

    Why: the builders used to write the UNIQUE dates next to ALL the image IDs
    ("1 date | 6 IDs"). Downstream pairs dates and IDs by position, so with
    more than one series per date the IDs were attached to the wrong dates.
    """
    by_date = {}
    for d, i in zip(dates, ids):
        if d is None or (not isinstance(d, str) and pd.isna(d)):
            continue
        key = str(d.date()) if hasattr(d, "date") else str(d)
        lst = by_date.setdefault(key, [])
        if i is None or (not isinstance(i, str) and pd.isna(i)):
            continue
        i = str(i).strip()
        if i.endswith(".0") and i[:-2].isdigit():
            i = i[:-2]
        if i and i.lower() not in ("nan", "none") and i not in lst:
            lst.append(i)
    sess, rep, allids = [], [], []
    for key, lst in by_date.items():
        sess.append(key)
        if not lst:
            rep.append("NA"); allids.append("NA"); continue
        nums = [x for x in lst if x.isdigit()]
        rep.append(min(nums, key=int) if len(nums) == len(lst) else lst[0])
        allids.append("/".join(lst))
    return sess, rep, allids


def session_id_columns(row: dict, T: str, dates, ids) -> None:
    """Write the aligned ID block of type T into `row` (dict of one subject):
    all_{T}_date, nb_total_{T}, {T}_IDs, all_{T}_image_ids, 1st_/Lst_{T}_ID."""
    sess, rep, allids = session_align_ids(dates, ids)
    real = [x for x in rep if x != "NA"]
    row[f"all_{T}_date"]      = "|".join(sess) if sess else np.nan
    row[f"nb_total_{T}"]      = len(sess)
    row[f"{T}_IDs"]           = "|".join(rep) if real else np.nan
    row[f"all_{T}_image_ids"] = "|".join(allids) if real else np.nan
    row[f"1st_{T}_ID"]        = real[0]  if real else np.nan
    row[f"Lst_{T}_ID"]        = real[-1] if real else np.nan
