#!/usr/bin/env python3
"""Score multi-species fruit-tree UniProt proteomes with PlantLLPS (grape-trained v13).

Reads a pinned species manifest (YAML/JSON), downloads FASTA in parallel (bounded),
scores each proteome with composition-proxy disorder, and deposits slim SQLite +
summary JSON (+ optional TSV.gz) under data/fruit_crops/.

Honesty: weights are grape-trained PlantLLPS v13 — not species-retrained.
Tier cutoffs are recalibrated on each proteome's score distribution.

Example::

    python3 scripts/score_fruit_crop_proteomes.py --list
    python3 scripts/score_fruit_crop_proteomes.py --species apple,peach
    python3 scripts/score_fruit_crop_proteomes.py --all
    python3 scripts/score_fruit_crop_proteomes.py --all --limit 500   # smoke
    python3 scripts/score_fruit_crop_proteomes.py --download-only --all
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sqlite3
import sys
import time
import urllib.parse
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import DEFAULT_UA  # noqa: E402
from plantllps.features import FEATURE_NAMES  # noqa: E402
from plantllps.model import PlantLLPSModel, assign_tier, percentile_cutoffs  # noqa: E402

DATA = ROOT / "data"
DEFAULT_MANIFEST = DATA / "fruit_crops" / "species_manifest.yaml"
DEFAULT_MODEL = DATA / "model_formal.json"
DEFAULT_OUT = DATA / "fruit_crops"

TSV_COLUMNS = [
    "protein_id",
    "length",
    "score_plantllps",
    "score_companion",
    "score_logistic",
    "tier",
    "percentile",
    *FEATURE_NAMES,
]


def load_manifest(path: Path) -> dict:
    text = path.read_text(encoding="utf-8")
    if path.suffix.lower() in {".yaml", ".yml"}:
        import yaml

        return yaml.safe_load(text)
    return json.loads(text)


def download_fasta(proteome: str, dest: Path, force: bool = False) -> dict:
    """Download UniProt FASTA.gz for a proteome. Returns status dict."""
    if dest.exists() and not force and dest.stat().st_size > 1000:
        return {
            "proteome": proteome,
            "path": str(dest),
            "bytes": dest.stat().st_size,
            "cached": True,
            "ok": True,
        }
    dest.parent.mkdir(parents=True, exist_ok=True)
    q = urllib.parse.urlencode(
        {
            "compressed": "true",
            "format": "fasta",
            "query": f"(proteome:{proteome})",
        }
    )
    url = f"https://rest.uniprot.org/uniprotkb/stream?{q}"
    req = urllib.request.Request(url, headers={"User-Agent": DEFAULT_UA})
    last_err: Exception | None = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(req, timeout=900) as resp, open(dest, "wb") as out:
                while True:
                    chunk = resp.read(1024 * 256)
                    if not chunk:
                        break
                    out.write(chunk)
            if dest.stat().st_size < 1000:
                raise RuntimeError(f"FASTA too small ({dest.stat().st_size} bytes)")
            return {
                "proteome": proteome,
                "path": str(dest),
                "bytes": dest.stat().st_size,
                "cached": False,
                "ok": True,
            }
        except Exception as e:
            last_err = e
            wait = 4 * (2**attempt)
            print(f"  [{proteome}] download failed ({e}); retry in {wait}s", flush=True)
            time.sleep(wait)
            if dest.exists() and dest.stat().st_size < 1000:
                dest.unlink(missing_ok=True)
    return {
        "proteome": proteome,
        "path": str(dest),
        "ok": False,
        "error": str(last_err),
    }


def open_fasta(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def parse_uniprot_fasta(path: Path, limit: int = 0) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    pid = None
    chunks: list[str] = []
    with open_fasta(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pid and chunks:
                    records.append((pid, "".join(chunks)))
                    if limit and len(records) >= limit:
                        return records
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 else header
                chunks = []
            else:
                chunks.append(line.strip())
        if pid and chunks and (not limit or len(records) < limit):
            records.append((pid, "".join(chunks)))
    return records


def write_slim_sqlite(
    db_path: Path,
    rows: list[dict],
    *,
    crop: dict,
    model_path: Path,
    honesty: str,
) -> None:
    db_path.parent.mkdir(parents=True, exist_ok=True)
    if db_path.exists():
        db_path.unlink()
    conn = sqlite3.connect(str(db_path))
    conn.executescript(
        """
        CREATE TABLE meta (key TEXT PRIMARY KEY, value TEXT NOT NULL);
        CREATE TABLE proteins (
            protein_id TEXT PRIMARY KEY,
            length INTEGER NOT NULL,
            score REAL NOT NULL,
            score_companion REAL,
            score_logistic REAL,
            tier TEXT NOT NULL,
            percentile REAL,
            gene_name TEXT,
            protein_name TEXT,
            grape_ortholog TEXT,
            grape_ortholog_gene TEXT,
            orthodb_id TEXT
        );
        CREATE INDEX idx_fc_score ON proteins(score DESC);
        CREATE INDEX idx_fc_tier ON proteins(tier);
        """
    )
    meta = {
        "crop": crop["slug"],
        "label": crop.get("label", crop["slug"]),
        "species": crop["species"],
        "proteome": crop["proteome"],
        "taxon": crop["taxon"],
        "n_scored": len(rows),
        "model": str(model_path.name),
        "disorder": "composition_proxy",
        "honesty": honesty,
        "layer": "fruit_crops",
    }
    for k, v in meta.items():
        conn.execute(
            "INSERT INTO meta(key, value) VALUES(?, ?)",
            (k, json.dumps(v) if not isinstance(v, str) else v),
        )
    conn.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?)",
        ("summary_json", json.dumps(meta)),
    )
    conn.executemany(
        """
        INSERT INTO proteins(
            protein_id, length, score, score_companion, score_logistic,
            tier, percentile
        ) VALUES (?, ?, ?, ?, ?, ?, ?)
        """,
        [
            (
                r["protein_id"],
                int(r["length"]),
                float(r["score_plantllps"]),
                float(r["score_companion"]),
                float(r["score_logistic"]),
                r["tier"],
                float(r["percentile"]) if r.get("percentile") != "" else None,
            )
            for r in rows
        ],
    )
    conn.commit()
    conn.close()


def score_species(
    crop: dict,
    *,
    out_dir: Path,
    model: PlantLLPSModel,
    model_path: Path,
    fasta: Path | None,
    force_download: bool,
    limit: int,
    min_length: int,
    tier_a_pct: float,
    tier_c_pct: float,
    write_tsv: bool,
) -> dict:
    slug = crop["slug"]
    proteome = crop["proteome"]
    fasta_path = fasta or (out_dir / f"{proteome}.fasta.gz")
    if fasta is None:
        print(f"[{slug}] Downloading {proteome} …", flush=True)
        dl = download_fasta(proteome, fasta_path, force=force_download)
        if not dl.get("ok"):
            return {
                "slug": slug,
                "species": crop["species"],
                "proteome": proteome,
                "ok": False,
                "stage": "download",
                "error": dl.get("error", "download failed"),
            }
        print(f"[{slug}] FASTA {dl.get('bytes')} bytes (cached={dl.get('cached')})", flush=True)

    print(f"[{slug}] Parsing {fasta_path} …", flush=True)
    records = parse_uniprot_fasta(fasta_path, limit=limit or 0)
    clean: list[tuple[str, str]] = []
    for pid, seq in records:
        if len(seq) < min_length:
            continue
        aa = sum(1 for c in seq.upper() if c in "ACDEFGHIKLMNPQRSTVWY")
        if aa / max(1, len(seq)) >= 0.95:
            clean.append((pid, seq))
    records = clean
    if not records:
        return {
            "slug": slug,
            "species": crop["species"],
            "proteome": proteome,
            "ok": False,
            "stage": "parse",
            "error": "no sequences after length/AA filter",
        }

    print(f"[{slug}] Scoring {len(records)} proteins …", flush=True)
    t0 = time.time()
    rows: list[dict] = []
    for i, (pid, seq) in enumerate(records, 1):
        scored = model.score(seq, use_metapredict=False)
        rows.append({"protein_id": pid, "percentile": "", **scored})
        if i % 2500 == 0:
            elapsed = time.time() - t0
            rate = i / max(elapsed, 1e-6)
            eta = (len(records) - i) / max(rate, 1e-6)
            print(
                f"  [{slug}] {i}/{len(records)} ({rate:.0f}/s, ETA {eta/60:.1f} min)",
                flush=True,
            )

    scores = [r["score_plantllps"] for r in rows]
    cutoffs = percentile_cutoffs(scores, pct_a=tier_a_pct, pct_c=tier_c_pct)
    order = sorted(range(len(rows)), key=lambda i: scores[i])
    ranks = [0.0] * len(rows)
    i = 0
    while i < len(order):
        j = i
        while j < len(order) and scores[order[j]] == scores[order[i]]:
            j += 1
        avg = (i + 1 + j) / 2.0
        for k in range(i, j):
            ranks[order[k]] = avg
        i = j
    n = len(rows)
    tier_counts = {"A": 0, "B": 0, "C": 0, "D": 0}
    for idx, r in enumerate(rows):
        feats = {k: r[k] for k in FEATURE_NAMES if k in r}
        r["tier"] = assign_tier(r["score_plantllps"], feats, cutoffs)
        r["percentile"] = round(100.0 * ranks[idx] / n, 3)
        tier_counts[r["tier"]] = tier_counts.get(r["tier"], 0) + 1

    honesty = (
        f"Fruit-crop condensatome for {crop['species']} ({proteome}). "
        "Scores use the grape-trained PlantLLPS v13 hybrid with composition-proxy disorder "
        "(not metapredict). Tier cutoffs recalibrated on this proteome's score distribution. "
        "Not a species-specific retrained model."
    )
    if limit:
        honesty += f" SUBSET deposit: --limit {limit} (not full proteome)."
    if crop.get("note"):
        honesty += f" Note: {crop['note']}"

    tsv_name = f"{slug}_condensatome.proteome.tsv.gz"
    tsv_path = out_dir / tsv_name
    if write_tsv:
        with gzip.open(tsv_path, "wt", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=TSV_COLUMNS, extrasaction="ignore", delimiter="\t"
            )
            writer.writeheader()
            for row in rows:
                writer.writerow({k: row.get(k, "") for k in TSV_COLUMNS})

    db_path = out_dir / f"{slug}_condensatome.sqlite"
    write_slim_sqlite(db_path, rows, crop=crop, model_path=model_path, honesty=honesty)

    elapsed = time.time() - t0
    summary = {
        "slug": slug,
        "crop": slug,
        "label": crop.get("label", slug),
        "species": crop["species"],
        "proteome": proteome,
        "taxon": crop["taxon"],
        "n_fasta_parsed": len(records),
        "n_scored": len(rows),
        "limit": limit or None,
        "min_length": min_length,
        "tiers": tier_counts,
        "tier_cutoffs": cutoffs,
        "wall_seconds": round(elapsed, 1),
        "model": model_path.name,
        "disorder": "composition_proxy",
        "tsv": tsv_name if write_tsv else None,
        "sqlite": db_path.name,
        "honesty": honesty,
        "mean_score": round(sum(scores) / max(1, len(scores)), 6),
        "note": crop.get("note"),
        "ok": True,
    }
    summary_path = out_dir / f"{slug}_condensatome.summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"[{slug}] done n={len(rows)} tiers={tier_counts} → {db_path.name}", flush=True)
    return summary


def write_index(manifest: dict, out_dir: Path, results: list[dict]) -> Path:
    by_slug = {r["slug"]: r for r in results if r.get("slug")}
    entries = []
    for sp in manifest.get("species") or []:
        slug = sp["slug"]
        status = sp.get("status", "score")
        entry = {
            "slug": slug,
            "label": sp.get("label", slug),
            "species": sp.get("species"),
            "taxon": sp.get("taxon"),
            "proteome": sp.get("proteome"),
            "status": status,
            "odb_key": sp.get("odb_key"),
            "note": sp.get("note"),
        }
        if status == "link":
            entry["deposit"] = sp.get("deposit")
            entry["ok"] = True
            # pull summary if present next to deposit
            dep = Path(sp.get("deposit") or "")
            if not dep.is_absolute():
                dep = ROOT / dep
            side = dep.with_name(dep.name.replace(".sqlite", ".summary.json"))
            if "second_crop" in str(dep) and side.exists():
                try:
                    sm = json.loads(side.read_text())
                    entry["n_scored"] = sm.get("n_scored")
                    entry["tiers"] = sm.get("tiers")
                except Exception:
                    pass
            elif slug == "grape":
                # grape summary from proteome_summary if present
                gs = DATA / "proteome_summary.json"
                if gs.exists():
                    try:
                        sm = json.loads(gs.read_text())
                        entry["n_scored"] = sm.get("n_scored") or sm.get("n_proteins")
                        entry["tiers"] = sm.get("tiers") or sm.get("tier_counts")
                    except Exception:
                        pass
        else:
            r = by_slug.get(slug) or {}
            entry.update(
                {
                    "ok": bool(r.get("ok")),
                    "n_scored": r.get("n_scored"),
                    "tiers": r.get("tiers"),
                    "sqlite": r.get("sqlite"),
                    "error": r.get("error"),
                    "stage": r.get("stage"),
                    "mean_score": r.get("mean_score"),
                    "wall_seconds": r.get("wall_seconds"),
                }
            )
        entries.append(entry)

    index = {
        "version": manifest.get("version", 1),
        "model": manifest.get("model"),
        "model_note": (manifest.get("model_note") or "").strip(),
        "honesty": (manifest.get("honesty") or "").strip(),
        "generated_at_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "n_scored_species": sum(1 for e in entries if e.get("status") == "score" and e.get("ok")),
        "n_failed_species": sum(
            1 for e in entries if e.get("status") == "score" and e.get("ok") is False
        ),
        "unavailable": manifest.get("unavailable") or [],
        "species": entries,
    }
    path = out_dir / "index.json"
    path.write_text(json.dumps(index, indent=2) + "\n")
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", default=str(DEFAULT_MANIFEST))
    ap.add_argument("--model", default=str(DEFAULT_MODEL))
    ap.add_argument("--out-dir", default="")
    ap.add_argument("--species", default="", help="Comma-separated slugs (default: none; use --all)")
    ap.add_argument("--all", action="store_true", help="Score every status=score species")
    ap.add_argument("--list", action="store_true", help="Print manifest species and exit")
    ap.add_argument("--download-only", action="store_true")
    ap.add_argument("--download-workers", type=int, default=3)
    ap.add_argument("--force-download", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--min-length", type=int, default=80)
    ap.add_argument("--tier-a-pct", type=float, default=95.0)
    ap.add_argument("--tier-c-pct", type=float, default=50.0)
    ap.add_argument("--no-tsv", action="store_true", help="Skip writing proteome TSV.gz")
    ap.add_argument("--fasta-dir", default="", help="Optional directory of cached FASTA.gz")
    args = ap.parse_args()

    manifest_path = Path(args.manifest)
    if not manifest_path.exists():
        raise SystemExit(f"Manifest not found: {manifest_path}")
    manifest = load_manifest(manifest_path)
    species = list(manifest.get("species") or [])
    out_dir = Path(args.out_dir) if args.out_dir else ROOT / (manifest.get("out_dir") or DEFAULT_OUT)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.list:
        for sp in species:
            print(
                f"{sp['slug']:22} {sp.get('status','?'):6} {sp.get('proteome','—')}  "
                f"{sp.get('species')}  (tax {sp.get('taxon')})"
            )
        print("\nUnavailable:")
        for u in manifest.get("unavailable") or []:
            print(f"  {u['slug']:22} {u.get('requested')} — {u.get('reason')}")
        return

    score_targets = [s for s in species if s.get("status") == "score"]
    if args.all:
        selected = score_targets
    elif args.species:
        want = {x.strip() for x in args.species.split(",") if x.strip()}
        selected = [s for s in score_targets if s["slug"] in want]
        missing = want - {s["slug"] for s in selected}
        # allow explicit link species? no — only score
        if missing:
            # maybe user named a skip/link
            known = {s["slug"] for s in species}
            bad = missing - known
            linked = missing & {s["slug"] for s in species if s.get("status") != "score"}
            if bad:
                raise SystemExit(f"Unknown slug(s): {sorted(bad)}")
            if linked:
                print(f"Skipping non-score slug(s): {sorted(linked)}", flush=True)
    else:
        raise SystemExit("Specify --all or --species apple,peach,...")

    if not selected:
        raise SystemExit("No species selected to score")

    # Parallel downloads
    fasta_dir = Path(args.fasta_dir) if args.fasta_dir else out_dir
    fasta_dir.mkdir(parents=True, exist_ok=True)
    download_jobs = []
    for sp in selected:
        dest = fasta_dir / f"{sp['proteome']}.fasta.gz"
        download_jobs.append((sp, dest))

    print(
        f"Downloading {len(download_jobs)} proteomes with {args.download_workers} workers …",
        flush=True,
    )
    dl_results: dict[str, dict] = {}
    with ThreadPoolExecutor(max_workers=max(1, args.download_workers)) as ex:
        futs = {
            ex.submit(download_fasta, sp["proteome"], dest, args.force_download): sp["slug"]
            for sp, dest in download_jobs
        }
        for fut in as_completed(futs):
            slug = futs[fut]
            try:
                dl_results[slug] = fut.result()
            except Exception as e:
                dl_results[slug] = {"ok": False, "error": str(e), "proteome": slug}
            status = dl_results[slug]
            if status.get("ok"):
                print(
                    f"  ✓ {slug} {status.get('proteome')} "
                    f"{status.get('bytes')} bytes cached={status.get('cached')}",
                    flush=True,
                )
            else:
                print(f"  ✗ {slug} {status.get('error')}", flush=True)

    if args.download_only:
        failures = [
            {"slug": s["slug"], "proteome": s["proteome"], **dl_results.get(s["slug"], {})}
            for s in selected
            if not (dl_results.get(s["slug"]) or {}).get("ok")
        ]
        report = {
            "download_only": True,
            "ok": [s["slug"] for s in selected if (dl_results.get(s["slug"]) or {}).get("ok")],
            "failed": failures,
            "unavailable": manifest.get("unavailable") or [],
        }
        (out_dir / "download_report.json").write_text(json.dumps(report, indent=2) + "\n")
        print(json.dumps(report, indent=2), flush=True)
        return

    model_path = Path(args.model)
    model = PlantLLPSModel.load(str(model_path))
    results: list[dict] = []
    for sp in selected:
        slug = sp["slug"]
        dl = dl_results.get(slug) or {}
        if not dl.get("ok"):
            results.append(
                {
                    "slug": slug,
                    "species": sp["species"],
                    "proteome": sp["proteome"],
                    "ok": False,
                    "stage": "download",
                    "error": dl.get("error", "download failed"),
                }
            )
            continue
        fasta = Path(dl["path"])
        try:
            summary = score_species(
                sp,
                out_dir=out_dir,
                model=model,
                model_path=model_path,
                fasta=fasta,
                force_download=False,
                limit=args.limit,
                min_length=args.min_length,
                tier_a_pct=args.tier_a_pct,
                tier_c_pct=args.tier_c_pct,
                write_tsv=not args.no_tsv,
            )
            results.append(summary)
        except Exception as e:
            results.append(
                {
                    "slug": slug,
                    "species": sp["species"],
                    "proteome": sp["proteome"],
                    "ok": False,
                    "stage": "score",
                    "error": str(e),
                }
            )
            print(f"[{slug}] SCORE FAILED: {e}", flush=True)

    index_path = write_index(manifest, out_dir, results)
    failures_path = out_dir / "failures.json"
    failures = [r for r in results if not r.get("ok")]
    # merge manifest unavailable
    failures_doc = {
        "score_failures": failures,
        "unavailable_no_proteome": manifest.get("unavailable") or [],
    }
    failures_path.write_text(json.dumps(failures_doc, indent=2) + "\n")

    ok_n = sum(1 for r in results if r.get("ok"))
    print(
        f"\nScored OK: {ok_n}/{len(results)}; index → {index_path}; "
        f"failures → {failures_path}",
        flush=True,
    )
    for r in results:
        if r.get("ok"):
            t = r.get("tiers") or {}
            print(
                f"  {r['slug']:22} n={r.get('n_scored')}  "
                f"A/B/C/D={t.get('A')}/{t.get('B')}/{t.get('C')}/{t.get('D')}",
                flush=True,
            )
        else:
            print(f"  {r['slug']:22} FAIL [{r.get('stage')}] {r.get('error')}", flush=True)


if __name__ == "__main__":
    main()
