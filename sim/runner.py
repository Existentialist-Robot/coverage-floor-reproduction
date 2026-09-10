"""
Shared run harness: deterministic seed derivation and a process-pool map.

Determinism contract: every run's seed is an explicit integer derived from
(experiment tag, cell key, replicate) by a stable hash, so a cell can be re-run in
isolation and reproduce bit-for-bit, and paired arms can be given IDENTICAL seeds
by construction ("matched seeds").  Scored seeds are always < 900000; calibration
seeds are >= 900000 and never enter a scored aggregate.
"""

from __future__ import annotations

import hashlib
import json
import os
import time
from concurrent.futures import ProcessPoolExecutor
from typing import Any, Callable, Iterable, Sequence

SCORED_SEED_CEILING = 900_000
HERE = os.path.dirname(os.path.abspath(__file__))
AGG = os.path.join(HERE, "aggregates")


def seed_for(tag: str, rep: int) -> int:
    """Stable scored seed for a (cell tag, replicate) pair.  Paired arms pass the
    SAME tag so they receive the same seed -- that is the matched-seed pairing."""
    h = hashlib.sha256(f"{tag}|{rep}".encode()).digest()
    return int.from_bytes(h[:6], "big") % SCORED_SEED_CEILING


def pmap(fn: Callable[[Any], Any], jobs: Sequence[Any], workers: int | None = None,
         label: str = "") -> list[Any]:
    workers = workers or max(1, (os.cpu_count() or 2) - 2)
    t0 = time.time()
    if workers == 1:
        res = [fn(j) for j in jobs]
    else:
        chunk = max(1, len(jobs) // (workers * 8))
        with ProcessPoolExecutor(max_workers=workers) as ex:
            res = list(ex.map(fn, jobs, chunksize=chunk))
    print(f"  {label}: {len(jobs)} runs, {workers} workers, "
          f"{time.time() - t0:.1f}s", flush=True)
    return res


def write_agg(name: str, payload: dict[str, Any]) -> str:
    os.makedirs(AGG, exist_ok=True)
    path = os.path.join(AGG, name)
    payload = dict(payload)
    payload.setdefault("_generated_by", "papers/bc-complex-networks-2026/sim")
    with open(path, "w") as f:
        json.dump(payload, f, indent=1, sort_keys=False, default=_default)
    print(f"  wrote aggregates/{name}", flush=True)
    return path


def _default(o):
    try:
        import numpy as np
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, np.ndarray):
            return o.tolist()
        if isinstance(o, (np.bool_,)):
            return bool(o)
    except Exception:
        pass
    raise TypeError(f"not serialisable: {type(o)}")


RULE_ENV = "SIM_DETECTION_RULE"


def load_rule(name: str | None = None) -> dict[str, Any]:
    """Load the pre-declared detection rule from configs/.

    Which file is read may be overridden by the SIM_DETECTION_RULE environment
    variable (set by each runner's --rule flag before any worker is spawned).  A
    spawned worker inherits the parent's environment, so parent and children always
    score under the SAME rule file -- that is what lets the whole grid be re-scored
    under the recalibrated rule (spec section 9, amendment A-8) without editing the
    rule the first campaign was scored under.
    """
    name = name or os.environ.get(RULE_ENV, "detection_rule.json")
    with open(os.path.join(HERE, "configs", name)) as f:
        return json.load(f)
