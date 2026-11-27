"""Shared schema + codec helpers for per-residue tracks (disorder / pLDDT).

Sidecar SQLite keeps the main grape atlas under the GitHub ~100 MiB limit while
still shipping NAR-reproducible offline deposits.
"""

from __future__ import annotations

import json
import sqlite3
import struct
import time
from pathlib import Path
from typing import Iterable, Optional, Sequence

SCHEMA_VERSION = 1

CREATE_SQL = """
CREATE TABLE IF NOT EXISTS meta (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS residue_tracks (
    protein_id TEXT PRIMARY KEY,
    length INTEGER NOT NULL,
    disorder_u8 BLOB,
    plddt_u8 BLOB,
    disorder_source TEXT,
    plddt_source TEXT,
    disorder_mean REAL,
    plddt_mean REAL,
    updated_at TEXT DEFAULT (datetime('now'))
);
CREATE INDEX IF NOT EXISTS idx_tracks_len ON residue_tracks(length);
"""


def open_tracks_db(path: Path, *, writable: bool = True) -> sqlite3.Connection:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    conn = sqlite3.connect(str(path), timeout=300.0)
    if writable:
        conn.execute("PRAGMA journal_mode=WAL")
        conn.execute("PRAGMA busy_timeout=120000")
        conn.execute("PRAGMA synchronous=NORMAL")
        conn.executescript(CREATE_SQL)
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("schema_version", str(SCHEMA_VERSION)),
        )
        conn.commit()
    return conn


def quantize_unit_interval(values: Sequence[float]) -> bytes:
    """Map [0,1] floats → uint8 (0–255)."""
    out = bytearray(len(values))
    for i, v in enumerate(values):
        x = float(v)
        if x < 0.0:
            x = 0.0
        elif x > 1.0:
            x = 1.0
        out[i] = int(round(x * 255.0))
    return bytes(out)


def dequantize_unit_interval(blob: bytes | memoryview | None) -> list[float]:
    if not blob:
        return []
    return [b / 255.0 for b in bytes(blob)]


def quantize_plddt(values: Sequence[float]) -> bytes:
    """Map pLDDT [0,100] → uint8."""
    out = bytearray(len(values))
    for i, v in enumerate(values):
        x = float(v)
        if x < 0.0:
            x = 0.0
        elif x > 100.0:
            x = 100.0
        out[i] = int(round(x))
    return bytes(out)


def dequantize_plddt(blob: bytes | memoryview | None) -> list[float]:
    if not blob:
        return []
    return [float(b) for b in bytes(blob)]


def upsert_track(
    conn: sqlite3.Connection,
    *,
    protein_id: str,
    length: int,
    disorder_u8: Optional[bytes] = None,
    plddt_u8: Optional[bytes] = None,
    disorder_source: Optional[str] = None,
    plddt_source: Optional[str] = None,
    disorder_mean: Optional[float] = None,
    plddt_mean: Optional[float] = None,
) -> None:
    """Insert/merge track columns without clobbering the other modality.

    Concurrent disorder + pLDDT depositors can share one DB: each call only
    overwrites fields it provides (COALESCE keeps the prior blob/source/mean).
    """
    conn.execute(
        """
        INSERT INTO residue_tracks(
            protein_id, length, disorder_u8, plddt_u8,
            disorder_source, plddt_source, disorder_mean, plddt_mean, updated_at
        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, datetime('now'))
        ON CONFLICT(protein_id) DO UPDATE SET
            length = CASE
                WHEN excluded.length > 0 THEN excluded.length
                ELSE residue_tracks.length
            END,
            disorder_u8 = COALESCE(excluded.disorder_u8, residue_tracks.disorder_u8),
            plddt_u8 = COALESCE(excluded.plddt_u8, residue_tracks.plddt_u8),
            disorder_source = COALESCE(excluded.disorder_source, residue_tracks.disorder_source),
            plddt_source = COALESCE(excluded.plddt_source, residue_tracks.plddt_source),
            disorder_mean = COALESCE(excluded.disorder_mean, residue_tracks.disorder_mean),
            plddt_mean = COALESCE(excluded.plddt_mean, residue_tracks.plddt_mean),
            updated_at = excluded.updated_at
        """,
        (
            protein_id,
            int(length or 0),
            disorder_u8,
            plddt_u8,
            disorder_source,
            plddt_source,
            disorder_mean,
            plddt_mean,
        ),
    )


def set_meta(conn: sqlite3.Connection, key: str, value) -> None:
    if not isinstance(value, str):
        value = json.dumps(value, sort_keys=True)
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        (key, value),
    )


def coverage(conn: sqlite3.Connection) -> dict:
    row = conn.execute(
        """
        SELECT
          COUNT(*) AS n,
          SUM(disorder_u8 IS NOT NULL) AS n_disorder,
          SUM(plddt_u8 IS NOT NULL) AS n_plddt,
          SUM(disorder_u8 IS NOT NULL AND plddt_u8 IS NOT NULL) AS n_both
        FROM residue_tracks
        """
    ).fetchone()
    return {
        "n_proteins": int(row[0] or 0),
        "n_with_disorder": int(row[1] or 0),
        "n_with_plddt": int(row[2] or 0),
        "n_with_both": int(row[3] or 0),
        "schema_version": SCHEMA_VERSION,
        "updated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def pack_f32_le(values: Iterable[float]) -> bytes:
    vals = list(values)
    return struct.pack(f"<{len(vals)}f", *vals)
