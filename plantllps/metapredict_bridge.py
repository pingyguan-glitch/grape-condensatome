"""Optional metapredict disorder scores (not required for core PlantLLPS).

Install (recommended dedicated env with PyTorch CPU/GPU)::

    micromamba create -n meta311 -c conda-forge python=3.11 pip
    micromamba run -n meta311 pip install metapredict

Then either:

* run scoring inside that env so ``import metapredict`` works, or
* set ``METAPREDICT_PYTHON`` to that interpreter and use
  ``scripts/score_metapredict_subset.py`` / ``--use-metapredict``.

The v13 logistic weights still expect the composition ``frac_disorder_proxy``.
When ``use_metapredict=True``, metapredict fraction is stored as
``frac_disorder_metapredict`` (and optionally used to *replace* the proxy for
companion composite / tiering via ``replace_proxy=True``).
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


def metapredict_available() -> bool:
    try:
        import metapredict  # noqa: F401

        return True
    except Exception:
        return False


def _meta_clean(sequence: str) -> str:
    """Keep canonical amino acids only (metapredict rejects X/U/O/B/Z)."""
    return "".join(ch for ch in sequence.upper() if ch in "ACDEFGHIKLMNPQRSTVWY")


def disorder_fraction(sequence: str, *, threshold: float = 0.5) -> float:
    """Fraction of residues with metapredict disorder score ≥ threshold."""
    import metapredict as meta

    seq = _meta_clean(sequence)
    if not seq:
        return 0.0
    scores = meta.predict_disorder(seq)
    n = len(scores)
    if n == 0:
        return 0.0
    return float(sum(1 for x in scores if float(x) >= threshold) / n)


def percent_disorder(sequence: str) -> float:
    """metapredict.percent_disorder (0–100) normalized to 0–1 when possible."""
    import metapredict as meta

    seq = _meta_clean(sequence)
    if not seq:
        return 0.0
    try:
        pct = float(meta.percent_disorder(seq))
        if pct > 1.0:
            return pct / 100.0
        return pct
    except Exception:
        return disorder_fraction(seq)


def score_records(
    records: Sequence[Tuple[str, str]],
    *,
    threshold: float = 0.5,
) -> Dict[str, float]:
    """Return protein_id → frac_disorder_metapredict for an in-process import."""
    out: Dict[str, float] = {}
    for pid, seq in records:
        out[pid] = disorder_fraction(seq, threshold=threshold)
    return out


def score_records_via_subprocess(
    records: Sequence[Tuple[str, str]],
    *,
    python_exe: Optional[str] = None,
    threshold: float = 0.5,
) -> Dict[str, float]:
    """Score via ``METAPREDICT_PYTHON`` (or ``python_exe``) when local import fails."""
    exe = python_exe or os.environ.get("METAPREDICT_PYTHON") or sys.executable
    with tempfile.TemporaryDirectory(prefix="metapredict_") as td:
        td_path = Path(td)
        fa = td_path / "in.fasta"
        out_json = td_path / "out.json"
        with fa.open("w") as fh:
            for pid, seq in records:
                clean = "".join(ch for ch in seq.upper() if ch in "ACDEFGHIKLMNPQRSTVWY")
                fh.write(f">{pid}\n")
                for i in range(0, len(clean), 80):
                    fh.write(clean[i : i + 80] + "\n")
        helper = r"""
import json, sys
from pathlib import Path
import metapredict as meta
fa, out, thr = Path(sys.argv[1]), Path(sys.argv[2]), float(sys.argv[3])
recs = []
pid = None
buf = []
for line in fa.read_text().splitlines():
    if line.startswith(">"):
        if pid is not None:
            recs.append((pid, "".join(buf)))
        pid = line[1:].split()[0]
        buf = []
    else:
        buf.append(line.strip())
if pid is not None:
    recs.append((pid, "".join(buf)))
result = {}
for p, s in recs:
    if not s:
        result[p] = 0.0
        continue
    scores = meta.predict_disorder(s)
    result[p] = float(sum(1 for x in scores if float(x) >= thr) / max(1, len(scores)))
out.write_text(json.dumps(result))
"""
        script = td_path / "run_meta.py"
        script.write_text(helper)
        proc = subprocess.run(
            [exe, str(script), str(fa), str(out_json), str(threshold)],
            capture_output=True,
            text=True,
            check=False,
        )
        if proc.returncode != 0 or not out_json.exists():
            raise RuntimeError(
                "metapredict subprocess failed "
                f"(exe={exe}, code={proc.returncode}):\n{proc.stderr[-2000:]}"
            )
        return {k: float(v) for k, v in json.loads(out_json.read_text()).items()}


def resolve_metapredict_scores(
    records: Sequence[Tuple[str, str]],
    *,
    threshold: float = 0.5,
    python_exe: Optional[str] = None,
) -> Dict[str, float]:
    if metapredict_available() and not python_exe:
        return score_records(records, threshold=threshold)
    return score_records_via_subprocess(
        records, python_exe=python_exe, threshold=threshold
    )


def attach_metapredict(
    feats: Mapping[str, float],
    meta_frac: float,
    *,
    replace_proxy: bool = False,
) -> Dict[str, float]:
    out = dict(feats)
    out["frac_disorder_metapredict"] = float(meta_frac)
    if replace_proxy:
        out["frac_disorder_proxy"] = float(meta_frac)
    return out
