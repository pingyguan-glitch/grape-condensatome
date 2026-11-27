"""GrapeCondensatome SQLite schema and helpers (prototype)."""

from __future__ import annotations

import gzip
import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional

_STRESS_CACHE: Optional[Dict[str, Any]] = None
_STRESS_SIDECAR = (
    Path(__file__).resolve().parents[1]
    / "data"
    / "expression_atlas"
    / "stress_by_protein.json.gz"
)


def _normalize_segments(segs: Any) -> List[Dict[str, Any]]:
    if not isinstance(segs, list):
        return []
    out: List[Dict[str, Any]] = []
    for s in segs:
        if not isinstance(s, dict):
            continue
        start = s.get("start", s.get("s"))
        end = s.get("end", s.get("e"))
        length = s.get("length")
        if length is None and start is not None and end is not None:
            try:
                length = int(end) - int(start) + 1
            except (TypeError, ValueError):
                length = None
        out.append(
            {
                **s,
                "start": start,
                "end": end,
                "length": length,
                "source": s.get("source") or s.get("src") or "",
                "label": s.get("label", s.get("lab")),
            }
        )
    return out


def _stress_from_sidecar(protein_id: Optional[str]) -> Optional[Dict[str, Any]]:
    global _STRESS_CACHE
    if not protein_id:
        return None
    if _STRESS_CACHE is None:
        _STRESS_CACHE = {}
        if _STRESS_SIDECAR.exists():
            try:
                with gzip.open(_STRESS_SIDECAR, "rt", encoding="utf-8") as fh:
                    _STRESS_CACHE = json.load(fh)
            except (OSError, json.JSONDecodeError):
                _STRESS_CACHE = {}
    hit = _STRESS_CACHE.get(protein_id)
    return hit if isinstance(hit, dict) else None


SCHEMA_SQL = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);

CREATE TABLE IF NOT EXISTS proteins (
    protein_id TEXT PRIMARY KEY,
    sequence_hash TEXT NOT NULL,
    length INTEGER NOT NULL,
    score REAL NOT NULL,
    score_companion REAL,
    score_logistic REAL,
    tier TEXT NOT NULL,
    features_json TEXT NOT NULL,
    annotations TEXT,
    source TEXT,
    created_at TEXT DEFAULT (datetime('now'))
);

CREATE INDEX IF NOT EXISTS idx_proteins_score ON proteins(score DESC);
CREATE INDEX IF NOT EXISTS idx_proteins_tier ON proteins(tier);
"""


def sequence_hash(seq: str) -> str:
    return hashlib.sha256(seq.upper().encode("ascii", errors="ignore")).hexdigest()[:16]


def connect(db_path: str | Path) -> sqlite3.Connection:
    db_path = Path(db_path)
    db_path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(db_path))
    conn.row_factory = sqlite3.Row
    conn.executescript(SCHEMA_SQL)
    return conn


def set_meta(conn: sqlite3.Connection, key: str, value: str) -> None:
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (key, value),
    )


def upsert_protein(
    conn: sqlite3.Connection,
    *,
    protein_id: str,
    sequence: str,
    score_row: Mapping[str, Any],
    annotations: Optional[str] = None,
    source: Optional[str] = None,
) -> None:
    feature_keys = [
        k
        for k in score_row.keys()
        if k
        not in {
            "protein_id",
            "score_plantllps",
            "score_companion",
            "score_logistic",
            "tier",
        }
    ]
    features = {k: score_row[k] for k in feature_keys}
    conn.execute(
        """
        INSERT INTO proteins(
            protein_id, sequence_hash, length, score, score_companion,
            score_logistic, tier, features_json, annotations, source
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
        ON CONFLICT(protein_id) DO UPDATE SET
            sequence_hash=excluded.sequence_hash,
            length=excluded.length,
            score=excluded.score,
            score_companion=excluded.score_companion,
            score_logistic=excluded.score_logistic,
            tier=excluded.tier,
            features_json=excluded.features_json,
            annotations=excluded.annotations,
            source=excluded.source,
            created_at=datetime('now')
        """,
        (
            protein_id,
            sequence_hash(sequence),
            int(score_row.get("length", len(sequence))),
            float(score_row["score_plantllps"]),
            float(score_row.get("score_companion", 0.0)),
            float(score_row.get("score_logistic", 0.0)),
            str(score_row["tier"]),
            json.dumps(features, sort_keys=True),
            annotations,
            source,
        ),
    )


def bulk_upsert(
    conn: sqlite3.Connection,
    rows: Iterable[tuple[str, str, Mapping[str, Any], Optional[str]]],
    source: str = "plantllps-prototype",
) -> int:
    n = 0
    for protein_id, sequence, score_row, annotations in rows:
        upsert_protein(
            conn,
            protein_id=protein_id,
            sequence=sequence,
            score_row=score_row,
            annotations=annotations,
            source=source,
        )
        n += 1
    conn.commit()
    return n


def summarize(conn: sqlite3.Connection) -> Dict[str, Any]:
    cur = conn.execute("SELECT COUNT(*) AS n FROM proteins")
    n = cur.fetchone()["n"]
    tiers = {
        r["tier"]: r["c"]
        for r in conn.execute("SELECT tier, COUNT(*) AS c FROM proteins GROUP BY tier")
    }
    top = [
        dict(r)
        for r in conn.execute(
            "SELECT protein_id, score, tier FROM proteins ORDER BY score DESC LIMIT 5"
        )
    ]
    return {"n_proteins": n, "tiers": tiers, "top": top}


def get_meta(conn: sqlite3.Connection) -> Dict[str, str]:
    return {r["key"]: r["value"] for r in conn.execute("SELECT key, value FROM meta")}


def get_protein(conn: sqlite3.Connection, protein_id: str) -> Optional[Dict[str, Any]]:
    row = conn.execute(
        "SELECT * FROM proteins WHERE protein_id = ? COLLATE NOCASE",
        (protein_id.strip(),),
    ).fetchone()
    if not row:
        return None
    return _row_to_protein(row)


def search_proteins(
    conn: sqlite3.Connection,
    query: str,
    *,
    limit: int = 50,
    offset: int = 0,
) -> tuple[list[Dict[str, Any]], int]:
    q = query.strip()
    if not q:
        return [], 0
    like = f"%{q}%"
    # Exact ID first, then substring on id / annotations
    exact = conn.execute(
        "SELECT * FROM proteins WHERE protein_id = ? COLLATE NOCASE",
        (q,),
    ).fetchone()
    if exact:
        return [_row_to_protein(exact)], 1
    count = conn.execute(
        """
        SELECT COUNT(*) AS n FROM proteins
        WHERE protein_id LIKE ? COLLATE NOCASE
           OR IFNULL(annotations, '') LIKE ? COLLATE NOCASE
           OR IFNULL(gene_name, '') LIKE ? COLLATE NOCASE
           OR IFNULL(protein_name, '') LIKE ? COLLATE NOCASE
           OR IFNULL(go_terms, '') LIKE ? COLLATE NOCASE
           OR IFNULL(keywords, '') LIKE ? COLLATE NOCASE
        """,
        (like, like, like, like, like, like),
    ).fetchone()["n"]
    rows = conn.execute(
        """
        SELECT * FROM proteins
        WHERE protein_id LIKE ? COLLATE NOCASE
           OR IFNULL(annotations, '') LIKE ? COLLATE NOCASE
           OR IFNULL(gene_name, '') LIKE ? COLLATE NOCASE
           OR IFNULL(protein_name, '') LIKE ? COLLATE NOCASE
           OR IFNULL(go_terms, '') LIKE ? COLLATE NOCASE
           OR IFNULL(keywords, '') LIKE ? COLLATE NOCASE
        ORDER BY score DESC
        LIMIT ? OFFSET ?
        """,
        (like, like, like, like, like, like, limit, offset),
    ).fetchall()
    return [_row_to_protein(r) for r in rows], int(count)


def list_by_tier(
    conn: sqlite3.Connection,
    tier: str,
    *,
    limit: int = 50,
    offset: int = 0,
    min_score: Optional[float] = None,
) -> tuple[list[Dict[str, Any]], int]:
    tier = tier.strip().upper()
    if tier not in {"A", "B", "C", "D"}:
        return [], 0
    if min_score is None:
        count = conn.execute(
            "SELECT COUNT(*) AS n FROM proteins WHERE tier = ?", (tier,)
        ).fetchone()["n"]
        rows = conn.execute(
            """
            SELECT * FROM proteins WHERE tier = ?
            ORDER BY score DESC LIMIT ? OFFSET ?
            """,
            (tier, limit, offset),
        ).fetchall()
    else:
        count = conn.execute(
            "SELECT COUNT(*) AS n FROM proteins WHERE tier = ? AND score >= ?",
            (tier, min_score),
        ).fetchone()["n"]
        rows = conn.execute(
            """
            SELECT * FROM proteins WHERE tier = ? AND score >= ?
            ORDER BY score DESC LIMIT ? OFFSET ?
            """,
            (tier, min_score, limit, offset),
        ).fetchall()
    return [_row_to_protein(r) for r in rows], int(count)


def top_proteins(conn: sqlite3.Connection, limit: int = 20) -> list[Dict[str, Any]]:
    rows = conn.execute(
        "SELECT * FROM proteins ORDER BY score DESC LIMIT ?",
        (limit,),
    ).fetchall()
    return [_row_to_protein(r) for r in rows]


def _row_to_protein(row: sqlite3.Row) -> Dict[str, Any]:
    d = dict(row)
    feats: Dict[str, Any] = {}
    raw = d.get("features_json") or "{}"
    try:
        feats = json.loads(raw)
    except json.JSONDecodeError:
        feats = {}
    d["features"] = feats
    d["percentile"] = feats.get("percentile")
    d["map65_tail_like"] = feats.get("map65_tail_like", 0)
    d["frac_disorder_proxy"] = feats.get("frac_disorder_proxy")
    d["frac_disorder_metapredict"] = d.get("frac_disorder_metapredict")
    if d["frac_disorder_metapredict"] is None:
        d["frac_disorder_metapredict"] = feats.get("frac_disorder_metapredict")
    ann_obj: Dict[str, Any] = {}
    ann_raw = d.get("annotations") or ""
    if ann_raw.startswith("{"):
        try:
            parsed = json.loads(ann_raw)
            if isinstance(parsed, dict):
                ann_obj = parsed
        except json.JSONDecodeError:
            ann_obj = {}
    d["annotation_obj"] = ann_obj
    d["gene_name"] = d.get("gene_name") or ann_obj.get("gene_name") or ""
    d["protein_name"] = d.get("protein_name") or ann_obj.get("protein_name") or ""
    go = d.get("go_terms") or ""
    if go:
        d["go_terms_list"] = [x for x in str(go).split(";") if x]
    else:
        d["go_terms_list"] = list(ann_obj.get("go_terms") or [])
    kw = d.get("keywords") or ""
    if kw:
        d["keywords_list"] = [x for x in str(kw).split(";") if x]
    else:
        d["keywords_list"] = list(ann_obj.get("keywords") or [])
    expr_raw = d.get("expr_tpm_json") or ""
    if expr_raw:
        try:
            d["expr_tpm"] = json.loads(expr_raw)
        except json.JSONDecodeError:
            d["expr_tpm"] = None
    else:
        d["expr_tpm"] = None
    stress_raw = d.get("expr_stress_json") or ""
    if stress_raw:
        try:
            d["expr_stress"] = json.loads(stress_raw)
        except json.JSONDecodeError:
            d["expr_stress"] = None
    else:
        d["expr_stress"] = _stress_from_sidecar(d.get("protein_id"))
    # AlphaFold / ortholog / comparator enrichments (columns may be absent on old DBs)
    segs_raw = d.get("af_disorder_segments_json") or ""
    if segs_raw:
        try:
            d["af_disorder_segments"] = _normalize_segments(json.loads(segs_raw))
        except json.JSONDecodeError:
            d["af_disorder_segments"] = []
    else:
        d["af_disorder_segments"] = []
    for key, dest in (
        ("domains_json", "domains"),
        ("idr_segments_json", "idr_segments"),
    ):
        raw = d.get(key) or ""
        if raw:
            try:
                d[dest] = _normalize_segments(json.loads(raw) or [])
            except json.JSONDecodeError:
                d[dest] = []
        else:
            d[dest] = []
    if not d["af_disorder_segments"] and d.get("idr_segments"):
        d["af_disorder_segments"] = [
            {
                "start": s.get("start"),
                "end": s.get("end"),
                "length": s.get("length"),
                "source": s.get("source"),
            }
            for s in d["idr_segments"]
            if "af_plddt" in str(s.get("source") or "")
        ]
    ortho_raw = d.get("ortholog_json") or ""
    if ortho_raw:
        try:
            d["ortholog_obj"] = json.loads(ortho_raw)
        except json.JSONDecodeError:
            d["ortholog_obj"] = None
    else:
        d["ortholog_obj"] = None
    if d.get("plaac_nllr_approx") is None:
        d["plaac_nllr_approx"] = feats.get("plaac_nllr_approx")
    return d
