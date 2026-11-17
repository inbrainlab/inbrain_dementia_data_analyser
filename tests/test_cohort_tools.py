import itertools
import random

import numpy as np

from inbrain_dementia_data_analyser.cohorts.longitudinal import _match_visits, _tolerances_days
from inbrain_dementia_data_analyser.cohorts.multimodal import _best_set


def _brute(events, required, window, min_mod):
    best = None
    labels = list(events)
    for k in range(min_mod, len(labels) + 1):
        for combo in itertools.combinations(labels, k):
            if any(r not in combo for r in required):
                continue
            for pick in itertools.product(*[events[l] for l in combo]):
                ts = [p[0] for p in pick]
                if max(ts) - min(ts) <= window:
                    key = (-len(combo), max(ts) - min(ts))
                    best = key if best is None or key < best else best
    return best


def test_multimodal_search_equals_exhaustive_search():
    rng = random.Random(1)
    for _ in range(300):
        events = {lab: sorted((float(rng.randint(0, 2000)), None) for _ in range(rng.randint(0, 4)))
                  for lab in ("A", "B", "C", "D")}
        events = {k: v for k, v in events.items() if v}
        required = [l for l in ("A", "B") if l in events][:rng.randint(0, 2)]
        res = _best_set(events, required, 365.0, max(2, len(required)), "most_modalities")
        brute = _brute(events, required, 365.0, max(2, len(required)))
        if brute is None:
            assert res is None
        else:
            assert res is not None and (res[0][0], res[0][1]) == brute


def test_tolerances_never_overlap():
    assert _tolerances_days([0, 365, 730], 200, clamp=True) == [182, 182, 182]
    assert _tolerances_days([0, 365], 100, clamp=True) == [100, 100]


def test_visit_matching_targets():
    sessions = [{"t": float(t), "acq": str(t), "id": f"I{t}", "all_ids": [f"I{t}"], "n_series": 1,
                 "manufacturer": None, "model": None} for t in (0, 30, 370, 800)]
    matches, unused = _match_visits(sessions, [0, 12, 24], 182, True, "visit")
    assert [m["id"] for m in matches] == ["I0", "I370", "I800"]
    assert len(unused) == 1
