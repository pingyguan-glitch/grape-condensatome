#!/usr/bin/env python3
"""Cache InterPro entry metadata (+ optional protein hit locations) for the atlas.

Offline deposit for NAR reproducibility:

* ``data/tracks/interpro_entries.json.gz`` — accession → {name, type, source_db}
  for every InterPro / Pfam id referenced in grape ``domain_*_ids``.
* Optional ``--protein-locations`` refreshes ``domains_json`` using InterPro
  API locations (all member databases) for gold ∪ Tier A∪B (or ``--limit``).

Fruit crops: pass ``--db`` / ``--out`` pointing at a fruit SQLite; same JSON schema.

Examples::

    python3 scripts/deposit_interpro_entries.py
    python3 scripts/deposit_interpro_entries.py --protein-locations --tier-ab --workers 12
"""

from __future__ import annotations

import argparse
import gzip
import json
import sqlite3
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import urlopen_text  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
DEFAULT_OUT = DATA / "tracks" / "interpro_entries.json.gz"
ENTRY_API = "https://www.ebi.ac.uk/interpro/api/entry/{db}/{acc}/"
PROTEIN_API = (
    "https://www.ebi.ac.uk/interpro/api/entry/all/protein/uniprot/{pid}/?page_size=100"
)


def collect_ids(conn: sqlite3.Connection) -> tuple[set[str], set[str]]:
    ipr: set[str] = set()
    pfam: set[str] = set()
    for (raw,) in conn.execute(
        "SELECT domain_interpro_ids FROM proteins "
        "WHERE domain_interpro_ids IS NOT NULL AND domain_interpro_ids != ''"
    ):
        for x in str(raw).split(";"):
            x = x.strip()
            if x:
                ipr.add(x)
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    if "domain_pfam_ids" in cols:
        for (raw,) in conn.execute(
            "SELECT domain_pfam_ids FROM proteins "
            "WHERE domain_pfam_ids IS NOT NULL AND domain_pfam_ids != ''"
        ):
            for x in str(raw).split(";"):
                x = x.strip()
                if x:
                    pfam.add(x)
    return ipr, pfam


def fetch_entry(db: str, acc: str) -> dict | None:
    url = ENTRY_API.format(db=db, acc=acc)
    try:
        raw = urlopen_text(url, timeout=45, accept="application/json")
    except Exception as e:
        code = getattr(e, "code", None)
        if code in (204, 404):
            return None
        raise
    data = json.loads(raw)
    meta = data.get("metadata") or data
    return {
        "accession": meta.get("accession") or acc,
        "name": meta.get("name") or meta.get("description") or acc,
        "type": meta.get("type") or "",
        "source_database": meta.get("source_database") or db,
    }


def fetch_protein_locations(pid: str) -> list[dict]:
    url = PROTEIN_API.format(pid=pid)
    try:
        raw = urlopen_text(url, timeout=60, accept="application/json")
    except Exception as e:
        code = getattr(e, "code", None)
        if code in (204, 404):
            return []
        raise
    data = json.loads(raw)
    domains: list[dict] = []
    for entry in data.get("results") or []:
        meta = entry.get("metadata") or {}
        acc = meta.get("accession")
        name = meta.get("name") or acc
        source_db = meta.get("source_database") or ""
        for prot in entry.get("proteins") or []:
            for loc in prot.get("entry_protein_locations") or []:
                for frag in loc.get("fragments") or []:
                    start = int(frag["start"])
                    end = int(frag["end"])
                    domains.append(
                        {
                            "start": start,
                            "end": end,
                            "length": end - start + 1,
                            "label": name,
                            "source": f"interpro_{source_db or 'entry'}",
                            "db": source_db or "InterPro",
                            "accession": acc,
                        }
                    )
    domains.sort(key=lambda d: (d["start"], d["end"], d.get("accession") or ""))
    return domains


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--out", default=str(DEFAULT_OUT))
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--protein-locations", action="store_true")
    ap.add_argument("--tier-ab", action="store_true", help="Locations only for Tier A∪B + gold")
    ap.add_argument("--limit", type=int, default=0, help="Max proteins for --protein-locations")
    ap.add_argument("--skip-entries", action="store_true")
    args = ap.parse_args()

    conn = sqlite3.connect(str(args.db), timeout=120.0)
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    entries: dict[str, dict] = {}
    if out_path.exists() and not args.force and not args.skip_entries:
        with gzip.open(out_path, "rt", encoding="utf-8") as fh:
            try:
                entries = json.load(fh).get("entries") or {}
            except json.JSONDecodeError:
                entries = {}

    report: dict = {"db": str(args.db), "out": str(out_path)}

    if not args.skip_entries:
        ipr, pfam = collect_ids(conn)
        jobs = [("interpro", a) for a in sorted(ipr) if a not in entries]
        jobs += [("pfam", a) for a in sorted(pfam) if a not in entries]
        print(f"Entry metadata pending: {len(jobs)} (have {len(entries)})", flush=True)
        t0 = time.time()
        n_ok = n_fail = 0
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futs = {pool.submit(fetch_entry, db, acc): (db, acc) for db, acc in jobs}
            done = 0
            for fut in as_completed(futs):
                db, acc = futs[fut]
                done += 1
                try:
                    meta = fut.result()
                except Exception:
                    n_fail += 1
                    meta = None
                if meta:
                    entries[acc] = meta
                    n_ok += 1
                else:
                    n_fail += 1
                if done % 500 == 0 or done == len(jobs):
                    print(
                        f"  entries {done}/{len(jobs)} ok={n_ok} fail={n_fail}",
                        flush=True,
                    )
        payload = {
            "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
            "n_entries": len(entries),
            "sources": {
                "interpro_api": "https://www.ebi.ac.uk/interpro/api/",
                "license_note": "InterPro data: EMBL-EBI; see https://www.ebi.ac.uk/interpro/about/interpro/",
            },
            "entries": entries,
        }
        with gzip.open(out_path, "wt", encoding="utf-8") as fh:
            json.dump(payload, fh, separators=(",", ":"), sort_keys=True)
        report["entries"] = {
            "n_ok": n_ok,
            "n_fail": n_fail,
            "n_total": len(entries),
            "elapsed_sec": round(time.time() - t0, 1),
        }
        print(f"Wrote {out_path} ({len(entries)} entries)", flush=True)

    if args.protein_locations:
        cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
        if "domains_json" not in cols:
            raise SystemExit("atlas DB lacks domains_json")
        if args.tier_ab:
            q = "SELECT protein_id FROM proteins WHERE tier IN ('A','B') ORDER BY score DESC"
        else:
            q = "SELECT protein_id FROM proteins ORDER BY score DESC"
        pids = [r[0] for r in conn.execute(q)]
        # Always include curated gold when present
        gold = DATA / "curated" / "llps_experimental_gold.tsv"
        if gold.exists():
            import csv

            have = set(pids)
            for row in csv.DictReader(gold.open(), delimiter="\t"):
                pid = (row.get("protein_id") or "").strip()
                if pid and pid not in have:
                    pids.insert(0, pid)
                    have.add(pid)
        if args.limit:
            pids = pids[: args.limit]
        print(f"InterPro protein locations: {len(pids)}", flush=True)
        t0 = time.time()
        n_ok = n_fail = n_upd = 0
        with ThreadPoolExecutor(max_workers=max(1, args.workers)) as pool:
            futs = {pool.submit(fetch_protein_locations, pid): pid for pid in pids}
            done = 0
            for fut in as_completed(futs):
                pid = futs[fut]
                done += 1
                try:
                    domains = fut.result()
                except Exception:
                    n_fail += 1
                    domains = None
                if domains is None:
                    continue
                n_ok += 1
                if domains:
                    # Merge with existing UniProt FT domains (keep both; de-dupe by start-end-acc)
                    existing = []
                    row = conn.execute(
                        "SELECT domains_json FROM proteins WHERE protein_id = ?", (pid,)
                    ).fetchone()
                    if row and row[0]:
                        try:
                            existing = json.loads(row[0]) or []
                        except json.JSONDecodeError:
                            existing = []
                    key = lambda d: (
                        int(d.get("start") or 0),
                        int(d.get("end") or 0),
                        str(d.get("accession") or ""),
                        str(d.get("source") or ""),
                    )
                    merged = {key(d): d for d in existing}
                    for d in domains:
                        merged[key(d)] = d
                    out_doms = sorted(merged.values(), key=lambda d: (d["start"], d["end"]))
                    conn.execute(
                        "UPDATE proteins SET domains_json = ? WHERE protein_id = ?",
                        (json.dumps(out_doms, sort_keys=True), pid),
                    )
                    n_upd += 1
                if done % 200 == 0 or done == len(pids):
                    conn.commit()
                    print(
                        f"  locations {done}/{len(pids)} ok={n_ok} fail={n_fail} upd={n_upd}",
                        flush=True,
                    )
        conn.commit()
        loc_summary = {
            "n_requested": len(pids),
            "n_ok": n_ok,
            "n_fail": n_fail,
            "n_updated": n_upd,
            "elapsed_sec": round(time.time() - t0, 1),
            "scope": "tier_ab_gold" if args.tier_ab else "all",
        }
        report["protein_locations"] = loc_summary
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            ("interpro_protein_locations", json.dumps(loc_summary)),
        )
        conn.commit()

    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("interpro_entries", json.dumps(report, sort_keys=True)),
    )
    conn.commit()
    report_path = out_path.with_suffix(".report.json")
    # .json.gz → write sibling .report.json next to stem
    report_path = Path(str(out_path).replace(".json.gz", ".report.json"))
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(json.dumps(report, indent=2), flush=True)
    conn.close()


if __name__ == "__main__":
    main()
