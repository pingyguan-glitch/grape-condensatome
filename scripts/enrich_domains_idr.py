#!/usr/bin/env python3
"""Join UniProt/InterPro domain annotations + IDR segment tracks into SQLite.

Sources (honest):
  * UniProt FT Domain / Pfam + InterPro xrefs for UP000009183 (stream cache).
  * Optional InterPro API backfill for Pfam match coordinates when UniProt
    Domain FT is empty but a Pfam xref exists (checkpointed).
  * IDR intervals from (1) existing AlphaFold low-pLDDT segments and
    (2) UniProt REGION \"Disordered\" (MobiDB-lite). Per-residue metapredict
    segments are **not** deposited (only mean frac), so they are not invented.

Example::

    python3 scripts/enrich_domains_idr.py
    python3 scripts/enrich_domains_idr.py --skip-interpro-backfill
"""

from __future__ import annotations

import argparse
import gzip
import json
import re
import sqlite3
import sys
import time
import urllib.parse
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import urlopen_text  # noqa: E402

DATA = ROOT / "data"
DEFAULT_DB = DATA / "grape_condensatome.sqlite"
CACHE = DATA / "domains" / "uniprot_up000009183_domains.tsv.gz"
INTERPRO_CKPT = DATA / "domains" / "interpro_pfam_backfill.checkpoint.json"

DOMAIN_RE = re.compile(
    r"DOMAIN\s+(\d+)\.\.(\d+);\s*/note=\"([^\"]*)\"",
    re.IGNORECASE,
)
REGION_RE = re.compile(
    r"REGION\s+(\d+)\.\.(\d+);\s*/note=\"([^\"]*)\"",
    re.IGNORECASE,
)


def ensure_columns(conn: sqlite3.Connection) -> None:
    cols = {r[1] for r in conn.execute("PRAGMA table_info(proteins)")}
    for col, decl in (
        ("domains_json", "TEXT"),
        ("domain_pfam_ids", "TEXT"),
        ("domain_interpro_ids", "TEXT"),
        ("idr_segments_json", "TEXT"),
        ("idr_n_residues", "INTEGER"),
        ("idr_sources", "TEXT"),
    ):
        if col not in cols:
            conn.execute(f"ALTER TABLE proteins ADD COLUMN {col} {decl}")


def download_uniprot(cache: Path, force: bool = False) -> Path:
    if cache.exists() and not force and cache.stat().st_size > 1000:
        print(f"Using cached UniProt domains → {cache}", flush=True)
        return cache
    cache.parent.mkdir(parents=True, exist_ok=True)
    fields = "accession,ft_domain,ft_region,xref_pfam,xref_interpro"
    qs = urllib.parse.urlencode(
        {
            "compressed": "true",
            "format": "tsv",
            "query": "(proteome:UP000009183)",
            "fields": fields,
        }
    )
    url = "https://rest.uniprot.org/uniprotkb/stream?" + qs
    print(f"Downloading UniProt domain/region annotations…\n  {url}", flush=True)
    raw = __import__("urllib.request").request.urlopen(
        __import__("urllib.request").request.Request(
            url, headers={"User-Agent": "GrapeCondensatome/0.6"}
        ),
        timeout=600,
    )
    with open(cache, "wb") as fh:
        while True:
            chunk = raw.read(1024 * 256)
            if not chunk:
                break
            fh.write(chunk)
    print(f"Wrote {cache} ({cache.stat().st_size} bytes)", flush=True)
    return cache


def split_ids(raw: str) -> list[str]:
    if not raw or not raw.strip():
        return []
    return [x.strip() for x in raw.replace(",", ";").split(";") if x.strip()]


def parse_domains_ft(raw: str) -> list[dict]:
    out = []
    for m in DOMAIN_RE.finditer(raw or ""):
        start, end = int(m.group(1)), int(m.group(2))
        out.append(
            {
                "start": start,
                "end": end,
                "length": end - start + 1,
                "label": m.group(3),
                "source": "uniprot_ft_domain",
                "db": None,
                "accession": None,
            }
        )
    return out


def parse_disordered_regions(raw: str) -> list[dict]:
    out = []
    for m in REGION_RE.finditer(raw or ""):
        note = m.group(3)
        if "disordered" not in note.lower():
            continue
        start, end = int(m.group(1)), int(m.group(2))
        out.append(
            {
                "start": start,
                "end": end,
                "length": end - start + 1,
                "source": "uniprot_mobidb",
                "label": note,
            }
        )
    return out


def load_uniprot_table(path: Path) -> dict[str, dict]:
    opener = gzip.open if str(path).endswith(".gz") else open
    out: dict[str, dict] = {}
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {h: i for i, h in enumerate(header)}

        def col(*names: str) -> int | None:
            for n in names:
                if n in idx:
                    return idx[n]
            return None

        i_acc = col("Entry", "accession")
        i_dom = col("Domain [FT]", "ft_domain")
        i_reg = col("Region", "ft_region")
        i_pfam = col("Pfam", "xref_pfam")
        i_ipr = col("InterPro", "xref_interpro")
        if i_acc is None:
            raise SystemExit(f"unexpected header: {header}")
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= i_acc:
                continue
            acc = parts[i_acc].strip()
            if not acc:
                continue
            dom_raw = parts[i_dom] if i_dom is not None and i_dom < len(parts) else ""
            reg_raw = parts[i_reg] if i_reg is not None and i_reg < len(parts) else ""
            pfam = split_ids(parts[i_pfam] if i_pfam is not None and i_pfam < len(parts) else "")
            ipr = split_ids(parts[i_ipr] if i_ipr is not None and i_ipr < len(parts) else "")
            domains = parse_domains_ft(dom_raw)
            # Attach sole Pfam id when exactly one domain + one pfam (best-effort)
            if len(domains) == 1 and len(pfam) == 1:
                domains[0]["db"] = "Pfam"
                domains[0]["accession"] = pfam[0]
            out[acc] = {
                "domains": domains,
                "pfam_ids": pfam,
                "interpro_ids": ipr,
                "disordered": parse_disordered_regions(reg_raw),
            }
    return out


def fetch_interpro_pfam(pid: str) -> list[dict]:
    url = (
        "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/"
        f"{pid}/?page_size=50"
    )
    try:
        raw = urlopen_text(url, timeout=45, accept="application/json")
    except Exception as e:
        code = getattr(e, "code", None)
        if code == 204 or code == 404:
            return []
        raise
    data = json.loads(raw)
    domains: list[dict] = []
    for entry in data.get("results") or []:
        meta = entry.get("metadata") or {}
        acc = meta.get("accession")
        name = meta.get("name") or acc
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
                            "source": "interpro_pfam",
                            "db": "Pfam",
                            "accession": acc,
                        }
                    )
    return domains


def merge_idr(af_segs: list[dict] | None, uniprot_disordered: list[dict]) -> tuple[list[dict], int, list[str]]:
    segs: list[dict] = []
    sources: set[str] = set()
    for s in af_segs or []:
        segs.append(
            {
                "start": int(s["start"]),
                "end": int(s["end"]),
                "length": int(s.get("length") or (int(s["end"]) - int(s["start"]) + 1)),
                "source": "af_plddt_lt50",
            }
        )
        sources.add("af_plddt_lt50")
    for s in uniprot_disordered:
        segs.append(dict(s))
        sources.add("uniprot_mobidb")
    segs.sort(key=lambda x: (x["start"], x["end"], x.get("source") or ""))
    # Residue coverage union (approximate unique residue count)
    covered: set[int] = set()
    for s in segs:
        covered.update(range(int(s["start"]), int(s["end"]) + 1))
    return segs, len(covered), sorted(sources)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--cache", default=str(CACHE))
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--skip-interpro-backfill", action="store_true")
    ap.add_argument("--interpro-workers", type=int, default=8)
    ap.add_argument("--interpro-limit", type=int, default=0)
    ap.add_argument("--checkpoint", default=str(INTERPRO_CKPT))
    args = ap.parse_args()

    cache = download_uniprot(Path(args.cache), force=args.force_download)
    table = load_uniprot_table(cache)
    print(f"UniProt domain rows: {len(table)}", flush=True)

    conn = sqlite3.connect(str(args.db))
    ensure_columns(conn)
    db_ids = [r[0] for r in conn.execute("SELECT protein_id FROM proteins")]
    print(f"SQLite proteins: {len(db_ids)}", flush=True)

    # InterPro backfill for Pfam-without-Domain-FT
    ckpt_path = Path(args.checkpoint)
    backfill: dict[str, list] = {}
    if ckpt_path.exists():
        try:
            backfill = json.loads(ckpt_path.read_text()).get("by_id") or {}
        except json.JSONDecodeError:
            backfill = {}

    need = []
    if not args.skip_interpro_backfill:
        for pid in db_ids:
            row = table.get(pid)
            if not row:
                continue
            if row["domains"]:
                continue
            if not row["pfam_ids"]:
                continue
            if pid in backfill:
                continue
            need.append(pid)
        if args.interpro_limit:
            need = need[: args.interpro_limit]
        print(f"InterPro Pfam backfill pending: {len(need)} (have {len(backfill)})", flush=True)
        t0 = time.time()
        done = 0
        with ThreadPoolExecutor(max_workers=max(1, args.interpro_workers)) as pool:
            futs = {pool.submit(fetch_interpro_pfam, pid): pid for pid in need}
            for fut in as_completed(futs):
                pid = futs[fut]
                try:
                    backfill[pid] = fut.result()
                except Exception as e:
                    backfill[pid] = []
                    print(f"  warn {pid}: {e}", flush=True)
                done += 1
                if done % 200 == 0 or done == len(need):
                    ckpt_path.parent.mkdir(parents=True, exist_ok=True)
                    ckpt_path.write_text(
                        json.dumps(
                            {
                                "by_id": backfill,
                                "n": len(backfill),
                                "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                            }
                        )
                    )
                    rate = done / max(time.time() - t0, 1e-6)
                    print(f"  interpro {done}/{len(need)} ({rate:.1f}/s)", flush=True)
        ckpt_path.parent.mkdir(parents=True, exist_ok=True)
        ckpt_path.write_text(
            json.dumps(
                {
                    "by_id": backfill,
                    "n": len(backfill),
                    "updated_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                }
            )
        )

    n_with_dom = 0
    n_with_idr = 0
    n_dom_entries = 0
    n_pfam_ids = 0
    for i, pid in enumerate(db_ids, 1):
        row = table.get(pid) or {
            "domains": [],
            "pfam_ids": [],
            "interpro_ids": [],
            "disordered": [],
        }
        domains = list(row["domains"])
        if not domains and pid in backfill:
            domains = list(backfill[pid])
        af_raw = conn.execute(
            "SELECT af_disorder_segments_json FROM proteins WHERE protein_id = ?",
            (pid,),
        ).fetchone()
        af_segs = []
        if af_raw and af_raw[0]:
            try:
                af_segs = json.loads(af_raw[0]) or []
            except json.JSONDecodeError:
                af_segs = []
        idr_segs, idr_n, idr_sources = merge_idr(af_segs, row["disordered"])
        pfam_ids = row["pfam_ids"]
        ipr_ids = row["interpro_ids"]
        if domains:
            n_with_dom += 1
            n_dom_entries += len(domains)
        if pfam_ids:
            n_pfam_ids += 1
        if idr_segs:
            n_with_idr += 1
        conn.execute(
            """
            UPDATE proteins SET
                domains_json = ?,
                domain_pfam_ids = ?,
                domain_interpro_ids = ?,
                idr_segments_json = ?,
                idr_n_residues = ?,
                idr_sources = ?
            WHERE protein_id = ?
            """,
            (
                json.dumps(domains, sort_keys=True) if domains else None,
                ";".join(pfam_ids) if pfam_ids else None,
                ";".join(ipr_ids) if ipr_ids else None,
                json.dumps(idr_segs, sort_keys=True) if idr_segs else None,
                idr_n if idr_segs else None,
                ",".join(idr_sources) if idr_sources else None,
                pid,
            ),
        )
        if i % 5000 == 0:
            conn.commit()
            print(f"  joined {i}/{len(db_ids)}", flush=True)

    summary = {
        "n_proteins": len(db_ids),
        "n_with_domain_coords": n_with_dom,
        "n_domain_entries": n_dom_entries,
        "n_with_pfam_xref": n_pfam_ids,
        "n_with_idr_segments": n_with_idr,
        "uniprot_cache": str(cache.name),
        "interpro_backfill_n": len(backfill),
        "idr_sources": ["af_plddt_lt50", "uniprot_mobidb"],
        "note": (
            "Domain coordinates from UniProt FT Domain and InterPro Pfam backfill; "
            "Pfam/InterPro ID lists from UniProt xrefs. IDR from AF pLDDT<50 segments "
            "(where present) + UniProt REGION Disordered (MobiDB-lite). "
            "metapredict contributes mean frac only — no per-residue IDR track."
        ),
    }
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("domains_idr", json.dumps(summary)),
    )
    conn.commit()
    print(json.dumps(summary, indent=2))
    for row in conn.execute(
        "SELECT protein_id, domain_pfam_ids, idr_n_residues, idr_sources, substr(domains_json,1,120) "
        "FROM proteins WHERE protein_id IN ('F6HET4','D7SSD8')"
    ):
        print(" ", row)
    conn.close()


if __name__ == "__main__":
    main()
