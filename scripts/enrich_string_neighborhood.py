#!/usr/bin/env python3
"""STRING PPI neighborhood + condensatome enrichment for the grape atlas.

Downloads *Vitis vinifera* (NCBI taxid 29760) STRING v12 protein links / aliases,
maps string IDs → UniProt accessions in ``grape_condensatome.sqlite``, then:

1. Builds an undirected graph of **intra-atlas** edges with combined_score ≥ threshold
   (default 400 = STRING medium confidence).
2. For every mapped protein with ≥1 neighbor, scores Tier A/B enrichment of the
   1-hop neighborhood (fold vs proteome background + hypergeometric p-value).
3. Writes:
   - ``string_neighborhood`` — **seed rows only** in SQLite (stats + top-N JSON)
   - ``data/string/neighborhood_stats.tsv.gz`` — full-proteome enrichment stats
   - ``data/string/top_neighbors.tsv.gz`` — top-N neighbors for every protein
   - ``data/string/seed_edges.tsv.gz`` — undirected seed-adjacent edges
     (not duplicated into SQLite — keeps ``grape_condensatome.sqlite`` under
     GitHub's 100 MB limit)

Honesty
-------
- STRING combined_score mixes experiments, databases, co-expression, textmining,
  neighborhood, fusion, and co-occurrence — **not** LLPS-specific evidence.
- Species is native Vitis STRING (29760), not Ath-transfer. Mapping uses UniProt_AC
  (+ preferred name / OrderedLocusName fallbacks); a few atlas IDs remain unmapped.
- Enrichment is descriptive adjacency bias, not causal condensate membership.

Example::

    python3 scripts/enrich_string_neighborhood.py
    python3 scripts/enrich_string_neighborhood.py --min-score 700 --seeds A
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sqlite3
import sys
import time
import urllib.request
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import DEFAULT_UA  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
CACHE = DATA / "string"

TAXON = 29760
STRING_VERSION = "12.0"
DOWNLOAD_BASE = "https://stringdb-downloads.org/download"

FILES = {
    "aliases": f"{TAXON}.protein.aliases.v{STRING_VERSION}.txt.gz",
    "info": f"{TAXON}.protein.info.v{STRING_VERSION}.txt.gz",
    "links": f"{TAXON}.protein.links.v{STRING_VERSION}.txt.gz",
}


def download(name: str, dest: Path, force: bool = False) -> Path:
    if dest.exists() and not force and dest.stat().st_size > 1000:
        print(f"Using cache {dest} ({dest.stat().st_size} bytes)", flush=True)
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    # paths: protein.aliases.v12.0/<file>, etc.
    kind = {
        FILES["aliases"]: f"protein.aliases.v{STRING_VERSION}",
        FILES["info"]: f"protein.info.v{STRING_VERSION}",
        FILES["links"]: f"protein.links.v{STRING_VERSION}",
    }[dest.name]
    url = f"{DOWNLOAD_BASE}/{kind}/{dest.name}"
    print(f"Downloading {url} → {dest}", flush=True)
    req = urllib.request.Request(url, headers={"User-Agent": DEFAULT_UA})
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=600) as resp, open(dest, "wb") as out:
                while True:
                    chunk = resp.read(1024 * 256)
                    if not chunk:
                        break
                    out.write(chunk)
            print(f"  wrote {dest.stat().st_size} bytes", flush=True)
            return dest
        except Exception as e:
            wait = 4 * (2**attempt)
            print(f"  failed ({e}); retry in {wait}s", flush=True)
            time.sleep(wait)
    raise SystemExit(f"download failed: {url}")


def ensure_files(force: bool = False) -> dict[str, Path]:
    out = {}
    for key, fname in FILES.items():
        out[key] = download(fname, CACHE / fname, force=force)
    return out


def build_mapping(
    aliases_path: Path,
    info_path: Path,
    atlas_ids: set[str],
    gene_to_id: dict[str, str],
) -> tuple[dict[str, str], dict[str, str]]:
    """Return string_id → uniprot, uniprot → string_id."""
    string_to_uni: dict[str, str] = {}
    uni_to_string: dict[str, str] = {}

    def add(sid: str, uni: str) -> None:
        if uni not in atlas_ids:
            return
        string_to_uni.setdefault(sid, uni)
        uni_to_string.setdefault(uni, sid)

    with gzip.open(aliases_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sid, alias, source = parts[0], parts[1], parts[2]
            if source in ("UniProt_AC", "Uniprot"):
                add(sid, alias)

    with gzip.open(info_path, "rt", encoding="utf-8", errors="replace") as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            sid, pref = parts[0], parts[1]
            if sid in string_to_uni:
                continue
            if pref in atlas_ids:
                add(sid, pref)
            elif pref.upper() in gene_to_id:
                add(sid, gene_to_id[pref.upper()])

    with gzip.open(aliases_path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            sid, alias, source = parts[0], parts[1], parts[2]
            if sid in string_to_uni:
                continue
            if source == "UniProt_GN_OrderedLocusNames" and alias.upper() in gene_to_id:
                add(sid, gene_to_id[alias.upper()])

    return string_to_uni, uni_to_string


def load_edges(
    links_path: Path,
    string_to_uni: dict[str, str],
    min_score: int,
) -> dict[str, dict[str, int]]:
    """Undirected adjacency: uni → {neighbor_uni: combined_score}."""
    adj: dict[str, dict[str, int]] = defaultdict(dict)
    n_raw = n_keep = 0
    with gzip.open(links_path, "rt", encoding="utf-8", errors="replace") as fh:
        next(fh, None)
        for line in fh:
            n_raw += 1
            parts = line.split()
            if len(parts) < 3:
                continue
            a, b, sc_s = parts[0], parts[1], parts[2]
            sc = int(sc_s)
            if sc < min_score:
                continue
            ua = string_to_uni.get(a)
            ub = string_to_uni.get(b)
            if not ua or not ub or ua == ub:
                continue
            n_keep += 1
            prev = adj[ua].get(ub)
            if prev is None or sc > prev:
                adj[ua][ub] = sc
                adj[ub][ua] = sc
    undirected = sum(len(v) for v in adj.values()) // 2
    print(
        f"Links: raw={n_raw:,} mapped_pair_rows={n_keep:,} "
        f"undirected≥{min_score}={undirected:,} proteins_with_neighbors={len(adj):,}",
        flush=True,
    )
    return adj


def hypergeom_sf(k: int, N: int, K: int, n: int) -> float:
    """P(X >= k) for X ~ Hypergeometric(N, K, n), clamped to (0, 1]."""
    if n <= 0 or N <= 0:
        return 1.0
    k = max(0, min(k, n, K))
    # survival sum; use log-space for stability on large N
    # P(X=i) = C(K,i)*C(N-K,n-i)/C(N,n)
    try:
        log_denom = math.lgamma(N + 1) - math.lgamma(n + 1) - math.lgamma(N - n + 1)
    except ValueError:
        return 1.0
    total = 0.0
    lo = max(0, n - (N - K))
    for i in range(max(k, lo), min(n, K) + 1):
        log_p = (
            math.lgamma(K + 1)
            - math.lgamma(i + 1)
            - math.lgamma(K - i + 1)
            + math.lgamma(N - K + 1)
            - math.lgamma(n - i + 1)
            - math.lgamma(N - K - (n - i) + 1)
            - log_denom
        )
        total += math.exp(log_p)
    if total <= 0:
        return 1.0
    return min(1.0, total)


def ensure_schema(con: sqlite3.Connection) -> None:
    con.executescript(
        """
        CREATE TABLE IF NOT EXISTS string_neighborhood (
            protein_id TEXT PRIMARY KEY,
            string_id TEXT,
            n_neighbors INTEGER,
            n_neighbors_tier_a INTEGER,
            n_neighbors_tier_ab INTEGER,
            frac_tier_ab REAL,
            enrichment_fold REAL,
            hypergeom_pvalue REAL,
            mean_neighbor_score REAL,
            mean_neighbor_combined REAL,
            top_neighbors_json TEXT,
            is_seed INTEGER DEFAULT 0,
            FOREIGN KEY (protein_id) REFERENCES proteins(protein_id)
        );
        CREATE INDEX IF NOT EXISTS idx_string_nb_enrich
            ON string_neighborhood(enrichment_fold DESC);
        CREATE INDEX IF NOT EXISTS idx_string_nb_pval
            ON string_neighborhood(hypergeom_pvalue ASC);

        """
    )
    # Drop legacy edge table if present (edges live in seed_edges.tsv.gz).
    con.execute("DROP TABLE IF EXISTS string_edges")
    # optional convenience columns on proteins (ignore if already present)
    cols = {r[1] for r in con.execute("PRAGMA table_info(proteins)")}
    for col, typ in (
        ("string_id", "TEXT"),
        ("string_n_neighbors", "INTEGER"),
        ("string_enrichment_fold", "REAL"),
        ("string_hypergeom_pvalue", "REAL"),
    ):
        if col not in cols:
            con.execute(f"ALTER TABLE proteins ADD COLUMN {col} {typ}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--min-score", type=int, default=400, help="STRING combined_score threshold")
    ap.add_argument(
        "--seeds",
        choices=("A", "AB"),
        default="A",
        help="High-tier seed set for edge deposition (default: Tier A)",
    )
    ap.add_argument(
        "--top-n",
        type=int,
        default=15,
        help="Top neighbors per protein in top_neighbors.tsv.gz (+ seed JSON in SQLite)",
    )
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--skip-download", action="store_true", help="Require cached STRING files")
    args = ap.parse_args()

    if not args.db.exists():
        raise SystemExit(f"DB not found: {args.db}")

    CACHE.mkdir(parents=True, exist_ok=True)
    if args.skip_download:
        paths = {k: CACHE / v for k, v in FILES.items()}
        missing = [p for p in paths.values() if not p.exists()]
        if missing:
            raise SystemExit(f"Missing STRING caches: {missing}")
    else:
        paths = ensure_files(force=args.force_download)

    con = sqlite3.connect(str(args.db))
    con.row_factory = sqlite3.Row
    rows = list(con.execute("SELECT protein_id, gene_name, tier, score FROM proteins"))
    atlas_ids = {r["protein_id"] for r in rows}
    tiers = {r["protein_id"]: r["tier"] for r in rows}
    scores = {r["protein_id"]: float(r["score"]) for r in rows}
    gene_to_id: dict[str, str] = {}
    for r in rows:
        g = (r["gene_name"] or "").strip()
        if g:
            gene_to_id[g.split()[0].upper()] = r["protein_id"]

    print(f"Atlas proteins: {len(atlas_ids):,}", flush=True)
    string_to_uni, uni_to_string = build_mapping(
        paths["aliases"], paths["info"], atlas_ids, gene_to_id
    )
    print(
        f"Mapped STRING↔UniProt: {len(uni_to_string):,} / {len(atlas_ids):,} "
        f"({100.0 * len(uni_to_string) / len(atlas_ids):.1f}%)",
        flush=True,
    )

    adj = load_edges(paths["links"], string_to_uni, args.min_score)

    n_ab = sum(1 for t in tiers.values() if t in ("A", "B"))
    n_all = len(atlas_ids)
    bg_rate = n_ab / n_all if n_all else 0.0
    seed_tiers = {"A"} if args.seeds == "A" else {"A", "B"}
    seeds = {pid for pid, t in tiers.items() if t in seed_tiers}
    print(f"Background Tier A∪B rate={bg_rate:.4f} (n={n_ab}); seeds={args.seeds} n={len(seeds)}", flush=True)

    ensure_schema(con)
    con.execute("DELETE FROM string_neighborhood")
    # clear convenience cols
    con.execute(
        "UPDATE proteins SET string_id=NULL, string_n_neighbors=NULL, "
        "string_enrichment_fold=NULL, string_hypergeom_pvalue=NULL"
    )

    gene_by_id = {r["protein_id"]: (r["gene_name"] or "") for r in rows}
    nb_rows = []  # all proteins — for TSV + seed SQLite subset
    top_tsv_rows: list[tuple] = []
    for pid, neighbors in adj.items():
        n = len(neighbors)
        n_a = sum(1 for x in neighbors if tiers.get(x) == "A")
        n_ab_n = sum(1 for x in neighbors if tiers.get(x) in ("A", "B"))
        frac = n_ab_n / n if n else 0.0
        fold = (frac / bg_rate) if bg_rate > 0 else None
        pval = hypergeom_sf(n_ab_n, n_all, n_ab, n)
        mean_llps = sum(scores.get(x, 0.0) for x in neighbors) / n if n else None
        mean_comb = sum(neighbors.values()) / n if n else None
        top = sorted(neighbors.items(), key=lambda kv: (-kv[1], -scores.get(kv[0], 0.0)))[
            : args.top_n
        ]
        is_seed = 1 if pid in seeds else 0
        top_json = json.dumps(
            [
                {
                    "protein_id": nid,
                    "combined_score": int(sc),
                    "tier": tiers.get(nid),
                    "score": round(scores.get(nid, 0.0), 6),
                    "gene_name": gene_by_id.get(nid) or "",
                }
                for nid, sc in top
            ],
            separators=(",", ":"),
        )
        for rank, (nid, sc) in enumerate(top, start=1):
            top_tsv_rows.append(
                (
                    pid,
                    nid,
                    int(sc),
                    rank,
                    tiers.get(nid) or "",
                    round(scores.get(nid, 0.0), 6),
                    gene_by_id.get(nid) or "",
                )
            )
        nb_rows.append(
            (
                pid,
                uni_to_string.get(pid),
                n,
                n_a,
                n_ab_n,
                frac,
                fold,
                pval,
                mean_llps,
                mean_comb,
                top_json,
                is_seed,
            )
        )

    # SQLite: seed rows only (GitHub 100 MB limit on grape_condensatome.sqlite)
    seed_nb = [r for r in nb_rows if r[11] == 1]
    con.executemany(
        """
        INSERT INTO string_neighborhood (
            protein_id, string_id, n_neighbors, n_neighbors_tier_a, n_neighbors_tier_ab,
            frac_tier_ab, enrichment_fold, hypergeom_pvalue, mean_neighbor_score,
            mean_neighbor_combined, top_neighbors_json, is_seed
        ) VALUES (?,?,?,?,?,?,?,?,?,?,?,?)
        """,
        seed_nb,
    )

    # deposit edges: undirected canonical order, ≥1 seed endpoint → TSV only
    edge_rows = []
    seen = set()
    for a in seeds:
        for b, sc in adj.get(a, {}).items():
            lo, hi = (a, b) if a < b else (b, a)
            key = (lo, hi)
            if key in seen:
                continue
            seen.add(key)
            edge_rows.append((lo, hi, int(sc)))

    # convenience columns for seeds + mapped string_id for all mapped proteins
    con.executemany(
        """
        UPDATE proteins SET
            string_id = ?,
            string_n_neighbors = ?,
            string_enrichment_fold = ?,
            string_hypergeom_pvalue = ?
        WHERE protein_id = ?
        """,
        [(r[1], r[2], r[6], r[7], r[0]) for r in seed_nb],
    )
    con.executemany(
        "UPDATE proteins SET string_id = ? WHERE protein_id = ?",
        [(sid, uni) for uni, sid in uni_to_string.items()],
    )

    summary = {
        "string_version_files": STRING_VERSION,
        "string_api_note": "Downloads from stringdb-downloads.org v12.0; API may report 12.5",
        "taxon_id": TAXON,
        "species": "Vitis vinifera",
        "min_combined_score": args.min_score,
        "score_meaning": (
            "STRING combined_score (0–1000); 150 low / 400 medium / 700 high / 900 highest. "
            "Channels mixed: experiments, databases, coexpression, textmining, "
            "neighborhood, fusion, cooccurrence — not LLPS-specific."
        ),
        "atlas_proteins": n_all,
        "mapped_to_string": len(uni_to_string),
        "mapping_frac": round(len(uni_to_string) / n_all, 6) if n_all else 0,
        "unmapped": n_all - len(uni_to_string),
        "proteins_with_neighbors": len(adj),
        "undirected_edges_ge_threshold": sum(len(v) for v in adj.values()) // 2,
        "seed_set": args.seeds,
        "n_seeds": len(seeds),
        "n_seed_edges_deposited": len(edge_rows),
        "n_neighborhood_rows_tsv": len(nb_rows),
        "n_neighborhood_rows_sqlite_seeds": len(seed_nb),
        "background_tier_ab": n_ab,
        "background_tier_ab_rate": round(bg_rate, 6),
        "enrichment_definition": (
            "fold = (frac neighbors in Tier A∪B) / (proteome Tier A∪B rate); "
            "p = Hypergeometric survival P(X>=k | N,K,n)"
        ),
        "limits": [
            "Native Vitis STRING network; no Ath/rice PPI transfer in this pass.",
            "Combined score ≠ experimental PPI proof; textmining can link co-mentioned genes.",
            "SQLite stores Tier A seed neighborhood rows only; full-proteome stats and "
            "top-N neighbors are in data/string/*.tsv.gz (GitHub 100 MB SQLite limit).",
            "Gold LLPS overlay (PR #7) is optional; seeds here are score tiers only.",
        ],
        "anchor_F6HET4": None,
    }
    if "F6HET4" in adj:
        row = next(r for r in nb_rows if r[0] == "F6HET4")
        top_anchor = [
            {
                "protein_id": nid,
                "combined_score": int(sc),
                "tier": tiers.get(nid),
                "score": round(scores.get(nid, 0.0), 6),
                "gene_name": gene_by_id.get(nid) or "",
            }
            for nid, sc in sorted(
                adj["F6HET4"].items(), key=lambda kv: (-kv[1], -scores.get(kv[0], 0.0))
            )[:10]
        ]
        summary["anchor_F6HET4"] = {
            "string_id": row[1],
            "n_neighbors": row[2],
            "n_tier_ab": row[4],
            "frac_tier_ab": row[5],
            "enrichment_fold": row[6],
            "hypergeom_pvalue": row[7],
            "top": top_anchor,
        }

    # meta keys
    con.execute(
        "CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT NOT NULL)"
    )
    meta = {
        "string_version": STRING_VERSION,
        "string_taxon": str(TAXON),
        "string_min_score": str(args.min_score),
        "string_mapped": str(len(uni_to_string)),
        "string_seed_set": args.seeds,
        "string_seed_edges": str(len(edge_rows)),
        "string_neighborhood_n": str(len(nb_rows)),
        "string_neighborhood_sqlite_seeds": str(len(seed_nb)),
    }
    for k, v in meta.items():
        con.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            (k, v),
        )

    con.commit()

    # deposit artifacts (small)
    CACHE.mkdir(parents=True, exist_ok=True)
    summary_path = CACHE / "enrichment_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    edges_path = CACHE / "seed_edges.tsv.gz"
    with gzip.open(edges_path, "wt", encoding="utf-8") as fh:
        fh.write("protein_a\tprotein_b\tcombined_score\ttier_a\ttier_b\n")
        for a, b, sc in sorted(edge_rows, key=lambda x: -x[2]):
            fh.write(f"{a}\t{b}\t{sc}\t{tiers.get(a,'')}\t{tiers.get(b,'')}\n")

    top_path = CACHE / "top_neighbors.tsv.gz"
    with gzip.open(top_path, "wt", encoding="utf-8") as fh:
        fh.write(
            "protein_id\tneighbor_id\tcombined_score\trank\tneighbor_tier\t"
            "neighbor_score\tneighbor_gene\n"
        )
        for row in top_tsv_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")

    stats_path = CACHE / "neighborhood_stats.tsv.gz"
    with gzip.open(stats_path, "wt", encoding="utf-8") as fh:
        fh.write(
            "protein_id\tstring_id\tn_neighbors\tn_tier_a\tn_tier_ab\tfrac_tier_ab\t"
            "enrichment_fold\thypergeom_pvalue\tmean_neighbor_score\t"
            "mean_neighbor_combined\tis_seed\ttier\tscore\n"
        )
        for r in nb_rows:
            fh.write(
                "\t".join(
                    [
                        r[0],
                        r[1] or "",
                        str(r[2]),
                        str(r[3]),
                        str(r[4]),
                        f"{r[5]:.6f}" if r[5] is not None else "",
                        f"{r[6]:.6f}" if r[6] is not None else "",
                        f"{r[7]:.6g}" if r[7] is not None else "",
                        f"{r[8]:.6f}" if r[8] is not None else "",
                        f"{r[9]:.6f}" if r[9] is not None else "",
                        str(r[11]),
                        tiers.get(r[0], ""),
                        f"{scores.get(r[0], 0.0):.6f}",
                    ]
                )
                + "\n"
            )

    readme = CACHE / "README.md"
    readme.write_text(
        f"""# STRING neighborhood deposit (Vitis)

- **Source:** STRING v{STRING_VERSION} downloads for NCBI taxid **{TAXON}** (*Vitis vinifera*)
  from `stringdb-downloads.org`.
- **Threshold:** combined_score ≥ **{args.min_score}** (medium by default).
- **Mapping:** UniProt_AC (+ preferred name / OrderedLocusName fallbacks) → atlas IDs.
  Coverage: **{len(uni_to_string):,} / {n_all:,}**.
- **SQLite seeds:** Tier **{args.seeds}** → **{len(seed_nb):,}** rows in `string_neighborhood`
  (top-N JSON included).
- **Full proteome (TSV.gz):** `neighborhood_stats.tsv.gz` ({len(nb_rows):,} rows),
  `top_neighbors.tsv.gz`, `seed_edges.tsv.gz` ({len(edge_rows):,} edges).
- **Enrichment:** Tier A∪B fraction among 1-hop neighbors vs proteome background
  ({bg_rate:.4f}); hypergeometric p-value. See `enrichment_summary.json`.

## Honesty

Combined scores mix experimental, database, co-expression, and **textmining** channels —
they are **not** LLPS evidence. This pass uses **native grape** STRING, not Arabidopsis
transfer. Large raw `*.protein.links*.txt.gz` caches are gitignored; re-download with
`python3 scripts/enrich_string_neighborhood.py`.

SQLite seed-only storage keeps `grape_condensatome.sqlite` under GitHub's 100 MB limit;
the Node UI loads the TSV.gz sidecars for non-seed protein pages.
""",
        encoding="utf-8",
    )

    con.execute("VACUUM")
    con.close()
    print(json.dumps(summary, indent=2)[:2000], flush=True)
    print(f"Wrote {summary_path}", flush=True)
    print(f"Wrote {edges_path} ({edges_path.stat().st_size} bytes)", flush=True)
    print(f"Wrote {top_path} ({top_path.stat().st_size} bytes)", flush=True)
    print(f"Wrote {stats_path} ({stats_path.stat().st_size} bytes)", flush=True)
    print(f"SQLite updated: {args.db} ({args.db.stat().st_size} bytes)", flush=True)


if __name__ == "__main__":
    main()
