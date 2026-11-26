"""GrapeCondensatome SQLite schema and helpers (prototype)."""

from __future__ import annotations

import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

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
    return d
