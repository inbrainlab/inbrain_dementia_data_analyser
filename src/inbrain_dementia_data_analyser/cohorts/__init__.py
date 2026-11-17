"""Cohort tools working on the expanded master (``MASTER_expanded_240m.csv``).

* :func:`create_cross_sectional_cohort` - one row per subject: anchor pairing
  with one event kept per subject, or a multimodal set;
* :func:`create_anchor_cohort` - anchor pairing with every option exposed;
* :func:`create_multimodal_cohort` - 2, 3, 4 ... modalities inside a window;
* :func:`create_longitudinal_cohort` - repeated acquisitions of one modality;
* :func:`create_longitudinal_multimodal_cohort` - a visit schedule where each
  visit carries a complete set of measures.
"""
from __future__ import annotations

import os

from .anchor import create_anchor_cohort
from .longitudinal import create_longitudinal_cohort
from .longitudinal_multimodal import create_longitudinal_multimodal_cohort
from .multimodal import DEFAULT_MODALITIES, create_multimodal_cohort
from .outputs import CohortResult

__all__ = ["CohortResult", "DEFAULT_MODALITIES", "create_anchor_cohort", "create_cross_sectional_cohort",
           "create_longitudinal_cohort", "create_longitudinal_multimodal_cohort",
           "create_multimodal_cohort"]


def create_cross_sectional_cohort(input_path: str | os.PathLike, *, method: str = "anchor",
                                  **kwargs) -> CohortResult:
    """One row per subject.

    Parameters
    ----------
    input_path:
        ``MASTER_expanded_240m.csv``.
    method:
        ``"anchor"`` - each subject's anchor acquisition (``anchor``) with the
        closest ``required`` / ``optional`` measures within ``margin_months``;
        one event per subject (``anchor_events``: ``"first"`` by default, or
        ``"last"`` / ``"most_optional"``). See :func:`create_anchor_cohort`.

        ``"multimodal"`` - the best set of ``modalities`` acquired within
        ``window_months``. See :func:`create_multimodal_cohort`.
    **kwargs:
        Passed to the chosen tool.

    Returns
    -------
    CohortResult
    """
    if method == "anchor":
        kwargs.setdefault("anchor_events", "first")
        if kwargs["anchor_events"] == "all":
            raise ValueError("anchor_events='all' keeps several events per subject; "
                             "use create_anchor_cohort for that")
        return create_anchor_cohort(input_path, **kwargs)
    if method == "multimodal":
        return create_multimodal_cohort(input_path, **kwargs)
    raise ValueError(f"method must be 'anchor' or 'multimodal', got {method!r}")
