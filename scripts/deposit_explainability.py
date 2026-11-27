#!/usr/bin/env python3
"""Deposit PlantLLPS SHAP / grammar attributions for gold + Tier A proteins.

Writes a gzip sidecar (keeps SQLite under GitHub's 100 MiB limit) plus a pretty
F6HET4 example JSON for NAR+ figure / protein-page demos.

Usage::

    python3 scripts/deposit_explainability.py
    python3 scripts/deposit_explainability.py --ids F6HET4 --no-occlusion  # smoke
"""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from plantllps.explain import explain_sequence  # noqa: E402
from plantllps.model import PlantLLPSModel  # noqa: E402

DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_FASTA = ROOT / "data" / "UP000009183.fasta.gz"
DEFAULT_MODEL = ROOT / "data" / "model_formal.json"
DEFAULT_OUT_DIR = ROOT / "data" / "explainability"
ANCHOR_ID = "F6HET4"


def load_fasta_ids(path: Path, want: Set[str]) -> Dict[str, str]:
    import gzip as gz

    opener = gz.open if str(path).endswith(".gz") else open
    found: Dict[str, str] = {}
    pid: Optional[str] = None
    chunks: List[str] = []
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:  # type: ignore[arg-type]
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None and pid in want:
                    found[pid] = "".join(chunks)
                    if len(found) >= len(want):
                        return found
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 and parts[0] in {"sp", "tr"} else header
                chunks = []
            else:
                if pid in want:
                    chunks.append(line)
        if pid is not None and pid in want:
            found[pid] = "".join(chunks)
    return found


def select_targets(con: sqlite3.Connection, extra_ids: Iterable[str]) -> List[sqlite3.Row]:
    con.row_factory = sqlite3.Row
    rows = con.execute(
        """
        SELECT protein_id, tier, score, llps_gold, length,
               idr_segments_json, domains_json
        FROM proteins
        WHERE llps_gold > 0 OR tier = 'A'
        ORDER BY CASE WHEN protein_id = ? THEN 0 ELSE 1 END, score DESC
        """,
        (ANCHOR_ID,),
    ).fetchall()
    by_id = {r["protein_id"]: r for r in rows}
    for pid in extra_ids:
        if pid and pid not in by_id:
            r = con.execute(
                """
                SELECT protein_id, tier, score, llps_gold, length,
                       idr_segments_json, domains_json
                FROM proteins WHERE protein_id = ? COLLATE NOCASE
                """,
                (pid,),
            ).fetchone()
            if r:
                by_id[r["protein_id"]] = r
    # stable: anchor first, then score desc
    ordered = sorted(
        by_id.values(),
        key=lambda r: (0 if r["protein_id"] == ANCHOR_ID else 1, -float(r["score"] or 0)),
    )
    return ordered


def parse_json_list(raw: Optional[str]) -> List[Dict[str, Any]]:
    if not raw:
        return []
    try:
        obj = json.loads(raw)
        return obj if isinstance(obj, list) else []
    except json.JSONDecodeError:
        return []


def slim_payload(full: Dict[str, Any]) -> Dict[str, Any]:
    """Compact per-ID record for the gzip sidecar."""
    axes = (full.get("grammar") or {}).get("axes") or {}
    axis_compact = {
        k: {
            "hybrid": round(float(v.get("hybrid", 0.0)), 6),
            "logistic_shap_prob": round(float(v.get("logistic_shap_prob", 0.0)), 6),
            "companion": round(float(v.get("companion", 0.0)), 6),
        }
        for k, v in axes.items()
    }
    occ = full.get("occlusion") or {}
    segs = full.get("segment_attributions") or {}
    return {
        "protein_id": full.get("protein_id"),
        "tier": full.get("tier"),
        "llps_gold": full.get("llps_gold"),
        "length": full.get("length"),
        "score_plantllps": round(float(full.get("score_plantllps") or 0), 6),
        "score_companion": round(float(full.get("score_companion") or 0), 6),
        "score_logistic": round(float(full.get("score_logistic") or 0), 6),
        "grammar_axes": axis_compact,
        "group_hybrid": {
            k: round(float(v), 6)
            for k, v in ((full.get("grammar") or {}).get("group_hybrid") or {}).items()
        },
        "top_features": full.get("top_features") or [],
        "occlusion": {
            "window": occ.get("window"),
            "step": occ.get("step"),
            "base_score": occ.get("base_score"),
            "windows": occ.get("windows") or [],
            "top_positive": occ.get("top_positive") or [],
        },
        "segment_attributions": {
            "idr": (segs.get("idr") or [])[:8],
            "domain": (segs.get("domain") or [])[:8],
        },
        "grammar_windows": full.get("grammar_windows") or [],
        "method_summary": full.get("method_summary"),
    }


def write_readme(path: Path, summary: Dict[str, Any]) -> None:
    f6 = summary.get("anchor_F6HET4") or {}
    axes = f6.get("grammar_axes") or {}
    path.write_text(
        f"""# PlantLLPS explainability deposits

SHAP-equivalent **Ser-rich / charge-block / IDR-disorder** attributions for
GrapeCondensatome proteins in the curated LLPS gold overlay **or** score **Tier A**
(always including anchor **F6HET4**).

## Method

- **Logistic head:** exact linear SHAP in logit space (`φ_i = w_i · z_i`).
- **Companion composite:** transparent weighted term decomposition
  (disorder 0.55 / LC 0.25 / blockiness 0.20 + MAP65/Ser bonuses).
- **Hybrid grammar axes:** average of logistic probability-local SHAP (grouped)
  and companion terms mapped onto the same axes.
- **Residue windows:** Gly occlusion ΔPlantLLPS (SHAP-equivalent leave-one-window-out);
  windows tagged `ser_rich` / `charge_block` when composition warrants; joined to
  deposited IDR / domain intervals.

Not KernelSHAP / TreeSHAP — the logistic model admits an exact closed form.
Occlusion is the sequence-space equivalent for local drivers.

## Files

| Path | Notes |
|------|-------|
| `attributions.json.gz` | All gold ∪ Tier A payloads (compact) |
| `F6HET4.json` | Pretty NAR+ example for VvMAP65-1 |
| `summary.json` | Counts + F6HET4 axis snapshot |
| `README.md` | This file |

## Coverage

- Proteins deposited: **{summary.get("n_proteins")}**
- Gold: **{summary.get("n_gold")}** · Tier A: **{summary.get("n_tier_a")}** · Union: **{summary.get("n_union")}**
- With occlusion windows: **{summary.get("n_with_occlusion")}**
- Model: `{summary.get("model")}` · window={summary.get("window")} step={summary.get("step")}

## Anchor F6HET4 (VvMAP65-1)

- Tier **{f6.get("tier")}** · PlantLLPS ≈ **{f6.get("score_plantllps")}** · gold={f6.get("llps_gold")}
- Grammar hybrid axes: Ser-rich ≈ **{(axes.get("ser_rich") or {}).get("hybrid")}** ·
  charge-block ≈ **{(axes.get("charge_block") or {}).get("hybrid")}** ·
  IDR/disorder ≈ **{(axes.get("idr_disorder") or {}).get("hybrid")}**

Regenerate::

```bash
python3 scripts/deposit_explainability.py
```

UI: protein pages load `attributions.json.gz` at startup (gold ∪ Tier A); live
`/score` uses the JS SHAP port for feature/grammar bars without occlusion.
""",
        encoding="utf-8",
    )


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--fasta", default=str(DEFAULT_FASTA))
    ap.add_argument("--model", default=str(DEFAULT_MODEL))
    ap.add_argument("--out-dir", default=str(DEFAULT_OUT_DIR))
    ap.add_argument("--ids", nargs="*", default=[], help="Restrict to these IDs (debug)")
    ap.add_argument("--window", type=int, default=25)
    ap.add_argument("--step", type=int, default=12)
    ap.add_argument("--top-k-windows", type=int, default=20)
    ap.add_argument("--no-occlusion", action="store_true")
    ap.add_argument("--limit", type=int, default=0, help="Max proteins (0 = all targets)")
    args = ap.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    model = PlantLLPSModel.load(args.model)
    con = sqlite3.connect(args.db)

    if args.ids:
        targets = []
        con.row_factory = sqlite3.Row
        for pid in args.ids:
            r = con.execute(
                """
                SELECT protein_id, tier, score, llps_gold, length,
                       idr_segments_json, domains_json
                FROM proteins WHERE protein_id = ? COLLATE NOCASE
                """,
                (pid,),
            ).fetchone()
            if r:
                targets.append(r)
            else:
                print(f"warning: {pid} not in DB", file=sys.stderr)
    else:
        targets = select_targets(con, [ANCHOR_ID])
    if args.limit > 0:
        targets = targets[: args.limit]

    want = {r["protein_id"] for r in targets}
    print(f"Loading {len(want)} sequences from {args.fasta} …", flush=True)
    seqs = load_fasta_ids(Path(args.fasta), want)
    missing = sorted(want - set(seqs))
    if missing:
        print(f"warning: {len(missing)} IDs missing from FASTA (e.g. {missing[:5]})", flush=True)

    by_id: Dict[str, Any] = {}
    t0 = time.time()
    n_occ = 0
    for i, row in enumerate(targets, 1):
        pid = row["protein_id"]
        seq = seqs.get(pid)
        if not seq:
            continue
        full = explain_sequence(
            model,
            seq,
            protein_id=pid,
            tier=str(row["tier"] or ""),
            llps_gold=int(row["llps_gold"] or 0),
            idr_segments=parse_json_list(row["idr_segments_json"]),
            domains=parse_json_list(row["domains_json"]),
            window=args.window,
            step=args.step,
            top_k_windows=args.top_k_windows,
            do_occlusion=not args.no_occlusion,
        )
        compact = slim_payload(full)
        by_id[pid] = compact
        if compact.get("occlusion", {}).get("windows"):
            n_occ += 1
        if pid == ANCHOR_ID:
            (out_dir / "F6HET4.json").write_text(
                json.dumps(full, indent=2), encoding="utf-8"
            )
        if i % 50 == 0 or i == len(targets):
            elapsed = time.time() - t0
            print(
                f"  [{i}/{len(targets)}] deposited={len(by_id)} "
                f"elapsed={elapsed:.1f}s ({elapsed/max(1,i):.2f}s/prot)",
                flush=True,
            )

    gz_path = out_dir / "attributions.json.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as fh:
        json.dump(by_id, fh, separators=(",", ":"))

    n_gold = sum(1 for v in by_id.values() if int(v.get("llps_gold") or 0) > 0)
    n_a = sum(1 for v in by_id.values() if str(v.get("tier")) == "A")
    summary = {
        "n_proteins": len(by_id),
        "n_gold": n_gold,
        "n_tier_a": n_a,
        "n_union": len(by_id),
        "n_with_occlusion": n_occ,
        "model": str(Path(args.model).as_posix()),
        "window": args.window,
        "step": args.step,
        "selection": "llps_gold>0 OR tier=A (+ F6HET4)",
        "sidecar": str(gz_path.relative_to(ROOT)),
        "anchor_example": "F6HET4.json",
        "anchor_F6HET4": by_id.get(ANCHOR_ID),
        "honesty": (
            "Logistic SHAP is exact in logit space for the linear head; "
            "hybrid grammar axes mix logistic φ with companion decomposition; "
            "window occlusion is a SHAP-equivalent local explanation, not KernelSHAP. "
            "Sidecar kept out of SQLite to respect GitHub 100 MiB limit."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    write_readme(out_dir / "README.md", summary)

    # Optional meta pointer (skipped if DB is read-only / near GitHub 100 MiB limit).
    # Node loads the gzip sidecar directly; meta is informational only.
    try:
        con.execute(
            "CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
        )
        con.execute(
            "INSERT OR REPLACE INTO meta(key, value) VALUES (?, ?)",
            (
                "explainability",
                json.dumps(
                    {
                        "n_proteins": len(by_id),
                        "sidecar": "data/explainability/attributions.json.gz",
                        "anchor": "data/explainability/F6HET4.json",
                        "selection": summary["selection"],
                    }
                ),
            ),
        )
        con.commit()
    except sqlite3.Error as exc:
        print(f"warning: could not write SQLite meta.explainability ({exc})", flush=True)
    finally:
        con.close()

    print(
        f"Wrote {len(by_id)} attributions → {gz_path} "
        f"({gz_path.stat().st_size/1e6:.2f} MB); F6HET4.json + summary.json",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
