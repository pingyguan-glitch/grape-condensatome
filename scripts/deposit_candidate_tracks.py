#!/usr/bin/env python3
"""Deposit AFDB pLDDT (+ optional metapredict) for fruit high-score candidates (NAR A2).

Reads ``data/fruit_crops/cross_species/conservative_candidates.tsv.gz`` (and optional
species-specific peaks), collects representative UniProt IDs for the deposited
panel species, and writes
``data/tracks/crop_candidates_residue_tracks.sqlite`` (same schema as grape tracks).

Grape proteome-wide tracks stay in ``grape_residue_tracks.sqlite``. This script
does **not** invent new species — only accessions already present in the
cross-species deposit.

Examples::

    # pLDDT for all conservative-candidate representatives (critical path)
    python3 scripts/deposit_candidate_tracks.py --plddt --workers 20

    # Top OrthoDB-dense species only (~8 including grape orthologs)
    python3 scripts/deposit_candidate_tracks.py --plddt --disorder \\
        --species-limit 8 --workers 16

    # Also include panel-contrast species-specific peaks (optional)
    python3 scripts/deposit_candidate_tracks.py --plddt --include-peaks
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
import time
from collections import Counter, defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from deposit_residue_tracks import (  # noqa: E402
    deposit_disorder,
    deposit_plddt,
    fetch_prediction,
)
from tracks_lib import coverage, open_tracks_db, set_meta  # noqa: E402
from _http_util import urlopen_text  # noqa: E402
import sqlite3

DATA = ROOT / "data"
CROSS = DATA / "fruit_crops" / "cross_species"
DEFAULT_CAND = CROSS / "conservative_candidates.tsv.gz"
DEFAULT_PEAKS = CROSS / "species_specific_peaks.tsv.gz"
DEFAULT_TRACKS = DATA / "tracks" / "crop_candidates_residue_tracks.sqlite"
AF_CACHE = DATA / "alphafold" / "candidates"


def open_tsv_gz(path: Path):
    return gzip.open(path, "rt", encoding="utf-8", newline="")


def species_columns(fieldnames: list[str]) -> list[str]:
    out = []
    for f in fieldnames:
        if f.startswith("id__"):
            out.append(f[len("id__") :])
    return out


def collect_ids(
    path: Path,
    *,
    only_high: bool = True,
    species_allow: set[str] | None = None,
) -> dict[str, set[str]]:
    """Return species_slug → set(protein_id)."""
    by_sp: dict[str, set[str]] = defaultdict(set)
    with open_tsv_gz(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        slugs = species_columns(reader.fieldnames or [])
        for row in reader:
            for slug in slugs:
                if species_allow is not None and slug not in species_allow:
                    continue
                pid = (row.get(f"id__{slug}") or "").strip()
                if not pid:
                    continue
                if only_high:
                    tier = (row.get(f"tier__{slug}") or "").strip().upper()
                    if tier not in {"A", "B"}:
                        continue
                by_sp[slug].add(pid)
            # grape ortholog pointer (accession in grape atlas)
            g = (row.get("grape_ortholog") or "").strip()
            if g and (species_allow is None or "grape" in species_allow):
                by_sp["grape"].add(g)
    return by_sp


def pick_top_species(by_sp: dict[str, set[str]], limit: int) -> set[str]:
    """Prefer species with the most high-score candidate IDs; always keep grape if present."""
    counts = Counter({s: len(ids) for s, ids in by_sp.items() if s != "grape"})
    top = [s for s, _ in counts.most_common(max(0, limit - (1 if "grape" in by_sp else 0)))]
    allow = set(top)
    if "grape" in by_sp:
        allow.add("grape")
    return allow


def resolve_atlas_db(slug: str) -> Path | None:
    if slug == "grape":
        return DATA / "grape_condensatome.sqlite"
    if slug in {"soybean", "tomato"}:
        return DATA / "second_crop" / f"{slug}_condensatome.sqlite"
    p = DATA / "fruit_crops" / f"{slug}_condensatome.sqlite"
    return p if p.exists() else None


def fetch_uniprot_fasta(pids: list[str], dest: Path) -> dict[str, str]:
    """Download sequences for candidate IDs (batched UniProt stream)."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    seqs: dict[str, str] = {}
    # Resume from existing FASTA
    if dest.exists():
        pid = None
        buf: list[str] = []
        for line in dest.read_text().splitlines():
            if line.startswith(">"):
                if pid is not None:
                    seqs[pid] = "".join(buf)
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 else header
                buf = []
            else:
                buf.append(line.strip())
        if pid is not None:
            seqs[pid] = "".join(buf)
    missing = [p for p in pids if p not in seqs]
    print(f"UniProt FASTA: have {len(seqs)}, missing {len(missing)}", flush=True)
    batch = 100
    for i in range(0, len(missing), batch):
        chunk = missing[i : i + batch]
        query = " OR ".join(f"accession:{p}" for p in chunk)
        from urllib.parse import urlencode

        qs = urlencode(
            {
                "query": query,
                "format": "fasta",
                "compressed": "false",
            }
        )
        url = "https://rest.uniprot.org/uniprotkb/stream?" + qs
        try:
            raw = urlopen_text(
                url,
                timeout=120,
                accept="text/plain",
                ua="GrapeCondensatome/1.8 (candidate tracks; research)",
            )
        except Exception as e:
            print(f"  FASTA batch fail ({e})", flush=True)
            continue
        pid = None
        buf: list[str] = []
        for line in raw.splitlines():
            if line.startswith(">"):
                if pid is not None:
                    seqs[pid] = "".join(buf)
                header = line[1:].split()[0]
                parts = header.split("|")
                pid = parts[1] if len(parts) > 1 else header
                buf = []
            else:
                buf.append(line.strip())
        if pid is not None:
            seqs[pid] = "".join(buf)
        print(f"  FASTA {min(i + batch, len(missing))}/{len(missing)}", flush=True)
        time.sleep(0.2)
    # rewrite cache
    with dest.open("w") as fh:
        for pid, seq in sorted(seqs.items()):
            fh.write(f">sp|{pid}|\n")
            for j in range(0, len(seq), 80):
                fh.write(seq[j : j + 80] + "\n")
    return seqs


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", default=str(DEFAULT_CAND))
    ap.add_argument("--peaks", default=str(DEFAULT_PEAKS))
    ap.add_argument("--include-peaks", action="store_true")
    ap.add_argument("--tracks-db", default=str(DEFAULT_TRACKS))
    ap.add_argument("--plddt", action="store_true")
    ap.add_argument("--disorder", action="store_true")
    ap.add_argument("--workers", type=int, default=16)
    ap.add_argument("--force", action="store_true")
    ap.add_argument(
        "--species-limit",
        type=int,
        default=0,
        help="If >0, keep grape + top N−1 species by candidate count",
    )
    ap.add_argument(
        "--species",
        default="",
        help="Comma-separated slug allow-list (overrides --species-limit)",
    )
    ap.add_argument("--batch-size", type=int, default=50)
    ap.add_argument("--cache-dir", default=str(AF_CACHE))
    ap.add_argument(
        "--fasta-cache",
        default=str(DATA / "tracks" / "candidate_sequences.fasta"),
    )
    args = ap.parse_args()
    if not args.plddt and not args.disorder:
        ap.error("Specify --plddt and/or --disorder")

    by_sp = collect_ids(Path(args.candidates), only_high=True)
    if args.include_peaks and Path(args.peaks).exists():
        extra = collect_ids(Path(args.peaks), only_high=True)
        for s, ids in extra.items():
            by_sp[s].update(ids)

    if args.species.strip():
        allow = {s.strip() for s in args.species.split(",") if s.strip()}
    elif args.species_limit > 0:
        allow = pick_top_species(by_sp, args.species_limit)
    else:
        allow = set(by_sp)

    by_sp = {s: ids for s, ids in by_sp.items() if s in allow}
    all_ids: list[str] = sorted({pid for ids in by_sp.values() for pid in ids})
    print(
        json.dumps(
            {
                "n_species": len(by_sp),
                "species": {s: len(ids) for s, ids in sorted(by_sp.items())},
                "n_unique_ids": len(all_ids),
            },
            indent=2,
        ),
        flush=True,
    )

    tracks = open_tracks_db(Path(args.tracks_db), writable=True)
    # Synthetic atlas connection: minimal table listing protein_ids for deposit_plddt
    # Prefer grape atlas for AF model_id when available; else empty model map.
    grape = sqlite3.connect(str(DATA / "grape_condensatome.sqlite"), timeout=120.0)
    # Build a temp atlas with union of IDs (for deposit helpers)
    tmp_atlas = sqlite3.connect(":memory:")
    tmp_atlas.execute(
        "CREATE TABLE proteins(protein_id TEXT PRIMARY KEY, af_model_id TEXT, "
        "af_mean_plddt REAL, frac_disorder_metapredict REAL)"
    )
    tmp_atlas.execute("CREATE TABLE meta(key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    for pid in all_ids:
        row = grape.execute(
            "SELECT af_model_id, af_mean_plddt FROM proteins WHERE protein_id = ?",
            (pid,),
        ).fetchone()
        mid = row[0] if row else None
        mean = row[1] if row else None
        # Also try fruit DBs for existence (no af columns usually)
        tmp_atlas.execute(
            "INSERT OR IGNORE INTO proteins(protein_id, af_model_id, af_mean_plddt) VALUES (?,?,?)",
            (pid, mid, mean),
        )
    tmp_atlas.commit()

    report: dict = {
        "role": "NAR A2 crop high-score candidate tracks",
        "species": {s: len(ids) for s, ids in sorted(by_sp.items())},
        "n_unique_ids": len(all_ids),
        "tracks_db": str(args.tracks_db),
        "candidates": str(args.candidates),
        "include_peaks": bool(args.include_peaks),
    }

    if args.plddt:
        # For fruit IDs without af_mean, deposit_plddt still tries prediction API when force
        # Widen selection: mark unknown IDs as needing fetch by setting force path
        report["plddt"] = deposit_plddt(
            atlas=tmp_atlas,
            tracks=tracks,
            ids=all_ids,
            workers=args.workers,
            force=args.force,
            cache_dir=Path(args.cache_dir),
        )
        # IDs without af_mean were filtered — re-run with force for remaining
        have = {
            r[0]
            for r in tracks.execute(
                "SELECT protein_id FROM residue_tracks WHERE plddt_u8 IS NOT NULL"
            )
        }
        missing = [p for p in all_ids if p not in have]
        if missing:
            print(f"pLDDT second pass for {len(missing)} IDs without prior AF mean…", flush=True)
            # Seed af_mean_plddt NULL but ensure they are attempted: use force
            report["plddt_pass2"] = deposit_plddt(
                atlas=tmp_atlas,
                tracks=tracks,
                ids=missing,
                workers=args.workers,
                force=True,
                cache_dir=Path(args.cache_dir),
            )

    if args.disorder:
        seqs = fetch_uniprot_fasta(all_ids, Path(args.fasta_cache))
        fasta_path = Path(args.fasta_cache)
        report["disorder"] = deposit_disorder(
            atlas=tmp_atlas,
            tracks=tracks,
            ids=[p for p in all_ids if p in seqs],
            fasta=fasta_path,
            batch_size=args.batch_size,
            python_exe=None,
            force=args.force,
        )

    cov = coverage(tracks)
    report["coverage"] = cov
    set_meta(tracks, "coverage", cov)
    set_meta(tracks, "last_deposit", report)
    # Species→count in meta; full ID lists in sidecar JSON
    id_index = {s: sorted(ids) for s, ids in by_sp.items()}
    set_meta(tracks, "species_counts", {s: len(v) for s, v in id_index.items()})
    idx_path = Path(args.tracks_db).with_suffix(".ids.json")
    idx_path.write_text(
        json.dumps(
            {"species": {s: len(v) for s, v in id_index.items()}, "ids": id_index},
            indent=2,
        )
    )
    set_meta(tracks, "ids_index_path", str(idx_path))
    tracks.commit()
    cov_path = Path(args.tracks_db).with_suffix(".coverage.json")
    cov_path.write_text(json.dumps(report, indent=2, sort_keys=True))
    print(json.dumps(report, indent=2), flush=True)
    tracks.close()
    grape.close()
    tmp_atlas.close()


if __name__ == "__main__":
    main()
