#!/usr/bin/env python3
"""Crosswalk grape LLPS candidates to public RNA-binding resources (NAR layer A3).

Joins UniProt UP000009183 proteins against:

1. **RBP2GO 2.0 Ath** RBPs via OrthoDB ``ortholog_ath`` (RBP2GO has no *Vitis*).
2. **ATtRACT** Ath gene / motif rows (sparse) via the same Ath orthologs; deposits
   motif family notes for matched proteins.
3. **Native UniProt** grape evidence: GO:0003723 and Keyword RNA-binding.
4. **Sequence RG / RGG motif** scans on the grape FASTA (stdlib).

Deposits under ``data/rbp/`` (sidecars; SQLite meta pointer only — keeps the
atlas DB under GitHub's ~100 MiB soft limit). UI loads
``annotations_by_protein.json.gz``.

Example::

    python3 scripts/import_rbp_sources.py --rbp2go-csv /tmp/RBP_AT.csv --fetch-attract
    python3 scripts/enrich_rbp_crosswalk.py
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import json
import re
import sqlite3
import sys
import zipfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterable

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_DB = ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_FASTA = ROOT / "data" / "UP000009183.fasta.gz"
DEFAULT_OUT = ROOT / "data" / "rbp"
CACHE = DEFAULT_OUT / "cache"

# Motif call thresholds (documented honesty — sequence grammar, not CLIP proof)
RGG_MIN = 2
RG_MIN = 8


def iter_fasta(path: Path) -> Iterable[tuple[str, str]]:
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:  # type: ignore[arg-type]
        pid: str | None = None
        chunks: list[str] = []
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if pid is not None:
                    yield pid, "".join(chunks)
                header = line[1:].strip()
                pid = header.split()[0]
                if "|" in pid:
                    parts = header.split("|")
                    if len(parts) >= 2 and parts[0] in {"sp", "tr"}:
                        pid = parts[1]
                chunks = []
            else:
                chunks.append(line)
        if pid is not None:
            yield pid, "".join(chunks)


def skip_hash_comments(fh) -> None:
    while True:
        pos = fh.tell()
        line = fh.readline()
        if not line:
            return
        if not line.startswith("#"):
            fh.seek(pos)
            return


def load_rbp2go(path: Path) -> dict[str, dict]:
    out: dict[str, dict] = {}
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        skip_hash_comments(fh)
        for r in csv.DictReader(fh):
            uid = (r.get("Uniprot_ID") or "").strip().upper()
            if not uid:
                continue
            try:
                score = float(r.get("RBP2GO_Score") or 0)
            except ValueError:
                score = 0.0
            try:
                n_ds = int(float(r.get("Nb_Datasets") or 0))
            except ValueError:
                n_ds = 0
            out[uid] = {
                "uniprot_id": uid,
                "entry_name": (r.get("Entry_Name") or "").strip(),
                "gene_name": (r.get("Gene_Name") or "").strip(),
                "rbp2go_score": score,
                "n_datasets": n_ds,
                "rbp_status": (r.get("RBP_status") or "").strip(),
                "has_rbd": str(r.get("has_RBD") or "").upper() in {"TRUE", "1", "YES"},
                "rbd_list": (r.get("RBD_list") or "").strip()[:240],
                "species": (r.get("Species") or "").strip(),
            }
    return out


def load_attract(zip_path: Path) -> tuple[dict[str, list[dict]], dict]:
    """Return Ath gene_id → motif rows + organism summary."""
    ath: dict[str, list[dict]] = defaultdict(list)
    orgs: Counter[str] = Counter()
    with zipfile.ZipFile(zip_path) as zf:
        name = "ATtRACT_db.txt"
        raw = zf.read(name).decode("utf-8", errors="replace")
    reader = csv.DictReader(io.StringIO(raw), delimiter="\t")
    # Duplicate Experiment_description columns appear in the header; DictReader
    # keeps the last — fine for our fields of interest.
    n_rows = 0
    for r in reader:
        n_rows += 1
        org = (r.get("Organism") or "").strip()
        orgs[org] += 1
        if org != "Arabidopsis_thaliana":
            continue
        gid = (r.get("Gene_id") or r.get("Gene_name") or "").strip().upper()
        if not gid:
            continue
        ath[gid].append(
            {
                "gene_id": gid,
                "gene_name": (r.get("Gene_name") or "").strip(),
                "motif": (r.get("Motif") or "").strip(),
                "family": (r.get("Family") or "").strip(),
                "len": (r.get("Len") or "").strip(),
                "pubmed": (r.get("Pubmed") or "").strip(),
                "matrix_id": (r.get("Matrix_id") or "").strip(),
            }
        )
    summary = {
        "n_rows_total": n_rows,
        "n_organisms": len(orgs),
        "top_organisms": orgs.most_common(12),
        "n_ath_rows": sum(len(v) for v in ath.values()),
        "n_ath_genes": len(ath),
        "plant_gap": (
            "ATtRACT is animal-dominated; Arabidopsis contributes few motif rows "
            "and Vitis is absent. Motif notes on grape proteins are Ath-ortholog "
            "transfers, not grape CLIP evidence."
        ),
    }
    return ath, summary


def resolve_attract_agi_to_uniprot(
    agi_ids: Iterable[str], cache_path: Path
) -> dict[str, list[str]]:
    """Map ATtRACT Ath AGI / gene symbols → UniProt accessions (cached TSV)."""
    sys.path.insert(0, str(ROOT / "scripts"))
    from _http_util import urlopen_text  # noqa: WPS433

    cached: dict[str, list[str]] = {}
    if cache_path.exists():
        with cache_path.open(newline="", encoding="utf-8") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                cached[r["query"].upper()] = [
                    x for x in (r.get("uniprot_ids") or "").split(";") if x
                ]

    out: dict[str, list[str]] = {}
    rows_out: list[dict] = []
    for q in sorted({a.upper() for a in agi_ids if a}):
        if q in cached:
            out[q] = cached[q]
            rows_out.append({"query": q, "uniprot_ids": ";".join(cached[q])})
            continue
        # Prefer gene_exact for AGI; also try gene name
        ids: list[str] = []
        for term in (f"gene_exact:{q}", f"gene:{q}"):
            url = (
                "https://rest.uniprot.org/uniprotkb/search?"
                f"query={term}+AND+organism_id:3702&fields=accession&format=tsv&size=10"
            )
            try:
                text = urlopen_text(url, timeout=60)
            except Exception as e:
                print(f"  UniProt resolve failed for {q}: {e}", flush=True)
                continue
            lines = [ln for ln in text.splitlines() if ln and not ln.startswith("Entry")]
            for ln in lines:
                acc = ln.split("\t")[0].strip().upper()
                if acc and acc not in ids:
                    ids.append(acc)
            if ids:
                break
        out[q] = ids
        rows_out.append({"query": q, "uniprot_ids": ";".join(ids)})
        print(f"  ATtRACT {q} → {ids or '—'}", flush=True)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with cache_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=["query", "uniprot_ids"], delimiter="\t")
        w.writeheader()
        for r in rows_out:
            w.writerow(r)
    return out


def scan_rg_motifs(seq: str, *, rgg_min: int = RGG_MIN, rg_min: int = RG_MIN) -> dict:
    s = seq.upper().replace("U", "X")
    # Overlapping-friendly counts
    n_rgg = len(re.findall(r"(?=RGG)", s))
    n_rg = len(re.findall(r"(?=RG)", s))
    # Longest RGG-run (consecutive RGG units)
    runs = [len(m.group(0)) // 3 for m in re.finditer(r"(?:RGG)+", s)]
    max_rgg_run = max(runs) if runs else 0
    rg_call = bool(n_rgg >= rgg_min or n_rg >= rg_min or max_rgg_run >= 2)
    note_parts = []
    if n_rgg:
        note_parts.append(f"RGG×{n_rgg}")
    if max_rgg_run >= 2:
        note_parts.append(f"max_RGG_run={max_rgg_run}")
    if n_rg:
        note_parts.append(f"RG×{n_rg}")
    return {
        "n_rgg": n_rgg,
        "n_rg": n_rg,
        "max_rgg_run": max_rgg_run,
        "rg_motif_call": rg_call,
        "rg_motif_note": "; ".join(note_parts) if note_parts else "",
    }


def write_tsv_gz(path: Path, fieldnames: list[str], rows: Iterable[dict]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    n = 0
    with gzip.open(path, "wt", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore", delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)
            n += 1
    return n


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--db", type=Path, default=DEFAULT_DB)
    ap.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--rbp2go-csv", type=Path, default=None)
    ap.add_argument("--attract-zip", type=Path, default=None)
    ap.add_argument("--rgg-min", type=int, default=RGG_MIN)
    ap.add_argument("--rg-min", type=int, default=RG_MIN)
    args = ap.parse_args()

    rgg_min, rg_min = args.rgg_min, args.rg_min

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    cache = out_dir / "cache"

    rbp2go_path = args.rbp2go_csv or (cache / "RBP_AT.csv")
    attract_path = args.attract_zip or (cache / "ATtRACT.zip")
    if not rbp2go_path.exists():
        raise SystemExit(
            f"Missing RBP2GO Ath CSV: {rbp2go_path}. Run scripts/import_rbp_sources.py first."
        )
    if not attract_path.exists():
        raise SystemExit(
            f"Missing ATtRACT.zip: {attract_path}. Run scripts/import_rbp_sources.py --fetch-attract."
        )
    if not args.fasta.exists():
        raise SystemExit(f"Missing FASTA: {args.fasta}")

    rbp2go = load_rbp2go(rbp2go_path)
    attract_ath, attract_summary = load_attract(attract_path)

    # Map Ath AGI locus (AT1G…) → UniProt for ATtRACT joins via gene symbols on RBP2GO
    agi_to_uniprot: dict[str, str] = {}
    for uid, rec in rbp2go.items():
        # Gene_Name may be "LOS1" not AGI; Alias sometimes holds AGI — skip if absent
        g = rec["gene_name"].upper()
        if re.match(r"^AT[1-5MC]G\d+", g):
            agi_to_uniprot[g] = uid
    # Also parse Ath gene ids from ATtRACT that look like AGI and map via UniProt later
    # Build reverse: Ath UniProt → grape protein_ids from SQLite
    con = sqlite3.connect(args.db)
    proteins = list(
        con.execute(
            """SELECT protein_id, ortholog_ath, ortholog_ath_gene, tier, score,
                      gene_name, protein_name, go_terms, keywords, features_json
               FROM proteins"""
        )
    )
    ath_to_grape: dict[str, list[str]] = defaultdict(list)
    grape_rows: dict[str, dict] = {}
    for (
        pid,
        ath,
        ath_gene,
        tier,
        score,
        gene_name,
        protein_name,
        go_terms,
        keywords,
        features_json,
    ) in proteins:
        grape_rows[pid] = {
            "protein_id": pid,
            "ortholog_ath": (ath or "").strip().upper(),
            "ortholog_ath_gene": (ath_gene or "").strip(),
            "tier": tier,
            "score": float(score) if score is not None else None,
            "gene_name": gene_name or "",
            "protein_name": protein_name or "",
            "go_terms": go_terms or "",
            "keywords": keywords or "",
            "features_json": features_json or "{}",
        }
        if ath:
            ath_to_grape[str(ath).strip().upper()].append(pid)

    # ATtRACT Ath genes are AGI IDs — resolve AGI → UniProt via UniProt gene mapping
    # from ortholog_ath_gene when present, else leave unmatched.
    gene_to_ath_uniprot: dict[str, set[str]] = defaultdict(set)
    for pid, row in grape_rows.items():
        if row["ortholog_ath"] and row["ortholog_ath_gene"]:
            gene_to_ath_uniprot[row["ortholog_ath_gene"].upper()].add(row["ortholog_ath"])
    # Also index RBP2GO gene names
    for uid, rec in rbp2go.items():
        if rec["gene_name"]:
            gene_to_ath_uniprot[rec["gene_name"].upper()].add(uid)

    print("Resolving ATtRACT Ath AGI → UniProt…", flush=True)
    queries: set[str] = set(attract_ath.keys())
    for motifs in attract_ath.values():
        for m in motifs:
            if m.get("gene_name"):
                queries.add(m["gene_name"].upper())
    agi_map = resolve_attract_agi_to_uniprot(
        queries, out_dir / "attract_ath_uniprot_map.tsv"
    )

    attract_by_uniprot: dict[str, list[dict]] = defaultdict(list)
    for agi, motifs in attract_ath.items():
        uids: set[str] = set()
        if agi in agi_to_uniprot:
            uids.add(agi_to_uniprot[agi])
        for q in {agi, *[m.get("gene_name", "").upper() for m in motifs]}:
            if not q:
                continue
            uids.update(agi_map.get(q, []))
            uids |= gene_to_ath_uniprot.get(q, set())
        for uid in uids:
            attract_by_uniprot[uid].extend(motifs)

    # Sequence motif scan
    print("Scanning RG/RGG motifs…", flush=True)
    motif_by_pid: dict[str, dict] = {}
    for pid, seq in iter_fasta(args.fasta):
        if pid in grape_rows:
            motif_by_pid[pid] = scan_rg_motifs(seq, rgg_min=rgg_min, rg_min=rg_min)
    print(f"  scanned {len(motif_by_pid)} sequences in atlas", flush=True)

    # Build per-protein annotations
    ann: dict[str, dict] = {}
    for pid, row in grape_rows.items():
        sources: list[str] = []
        notes: list[str] = []
        rbp2go_hit = None
        attract_motifs: list[dict] = []
        ath = row["ortholog_ath"]
        if ath and ath in rbp2go:
            rbp2go_hit = rbp2go[ath]
            sources.append("rbp2go_ath")
            notes.append(
                f"RBP2GO Ath {rbp2go_hit['uniprot_id']} "
                f"({rbp2go_hit['gene_name'] or '—'}; score={rbp2go_hit['rbp2go_score']}; "
                f"datasets={rbp2go_hit['n_datasets']})"
            )
            if ath in attract_by_uniprot:
                attract_motifs = attract_by_uniprot[ath][:12]
            # Also try gene-name join for ATtRACT AGI
        # Gene-symbol path for ATtRACT when UniProt map missed
        if not attract_motifs and ath:
            gsym = (row["ortholog_ath_gene"] or "").upper()
            if gsym and gsym in attract_ath:
                attract_motifs = attract_ath[gsym][:12]
            # AGI-like gene names on Ath side via RBP2GO gene
            if rbp2go_hit and rbp2go_hit["gene_name"].upper() in attract_ath:
                attract_motifs = attract_ath[rbp2go_hit["gene_name"].upper()][:12]

        if attract_motifs:
            sources.append("attract_ath")
            fams = sorted({m["family"] for m in attract_motifs if m.get("family") and m["family"] != "N/A"})
            mots = sorted({m["motif"] for m in attract_motifs if m.get("motif")})[:6]
            notes.append(
                "ATtRACT Ath motifs: "
                + (", ".join(mots) if mots else "—")
                + (f" · families {', '.join(fams)}" if fams else "")
            )

        go = row["go_terms"]
        kw = row["keywords"]
        uniprot_go = "GO:0003723" in go
        uniprot_kw = "RNA-binding" in kw or "RNA binding" in kw
        if uniprot_go:
            sources.append("uniprot_go3723")
            notes.append("UniProt GO:0003723 RNA binding")
        if uniprot_kw:
            sources.append("uniprot_kw_rna_binding")
            notes.append("UniProt keyword RNA-binding")

        motif = motif_by_pid.get(pid) or {
            "n_rgg": 0,
            "n_rg": 0,
            "max_rgg_run": 0,
            "rg_motif_call": False,
            "rg_motif_note": "",
        }
        if motif["rg_motif_call"]:
            sources.append("rg_motif")
            notes.append(f"RG/RGG grammar ({motif['rg_motif_note']})")

        # RNA-related scoring context (not evidence by itself)
        try:
            feats = json.loads(row["features_json"] or "{}")
        except json.JSONDecodeError:
            feats = {}
        frac_arg = float(feats.get("frac_arg") or 0)
        high_arg = frac_arg >= 0.08
        tier_ab = row["tier"] in {"A", "B"}

        if not sources:
            continue

        ann[pid] = {
            "protein_id": pid,
            "gene_name": row["gene_name"],
            "protein_name": row["protein_name"],
            "tier": row["tier"],
            "score": row["score"],
            "ortholog_ath": ath or None,
            "ortholog_ath_gene": row["ortholog_ath_gene"] or None,
            "rbp_flag": 1,
            "sources": sources,
            "source_tags": ";".join(sources),
            "notes": notes,
            "motif_note": "; ".join(notes[:4]),
            "rg_motif": motif,
            "rbp2go": rbp2go_hit,
            "attract_motifs": attract_motifs,
            "uniprot_go_rna_binding": uniprot_go,
            "uniprot_kw_rna_binding": uniprot_kw,
            "frac_arg": frac_arg,
            "rna_score_context": bool(tier_ab and (high_arg or motif["rg_motif_call"])),
        }

    # Stats
    n_proteins = len(grape_rows)
    n_rbp = len(ann)
    by_src = Counter()
    for a in ann.values():
        for s in a["sources"]:
            by_src[s] += 1
    tier_c = Counter(a["tier"] for a in ann.values())
    n_rg_call = sum(1 for a in ann.values() if a["rg_motif"]["rg_motif_call"])
    n_rbp2go_grape = sum(1 for a in ann.values() if "rbp2go_ath" in a["sources"])
    n_attract_grape = sum(1 for a in ann.values() if "attract_ath" in a["sources"])
    n_tier_ab_rbp = sum(1 for a in ann.values() if a["tier"] in {"A", "B"})
    n_rg_and_rbp2go = sum(
        1
        for a in ann.values()
        if a["rg_motif"]["rg_motif_call"] and "rbp2go_ath" in a["sources"]
    )

    # Deposit slim RBP2GO table
    rbp2go_rows = [
        {
            "uniprot_id": r["uniprot_id"],
            "entry_name": r["entry_name"],
            "gene_name": r["gene_name"],
            "rbp2go_score": f"{r['rbp2go_score']:.4g}",
            "n_datasets": r["n_datasets"],
            "rbp_status": r["rbp_status"],
            "has_rbd": int(r["has_rbd"]),
            "rbd_list": r["rbd_list"],
            "species": r["species"],
            "n_grape_orthologs": len(ath_to_grape.get(r["uniprot_id"], [])),
        }
        for r in sorted(rbp2go.values(), key=lambda x: -x["rbp2go_score"])
    ]
    write_tsv_gz(
        out_dir / "rbp2go_arabidopsis.tsv.gz",
        [
            "uniprot_id",
            "entry_name",
            "gene_name",
            "rbp2go_score",
            "n_datasets",
            "rbp_status",
            "has_rbd",
            "rbd_list",
            "species",
            "n_grape_orthologs",
        ],
        rbp2go_rows,
    )

    # ATtRACT Ath motif deposit
    attract_rows = []
    for gid, motifs in sorted(attract_ath.items()):
        for m in motifs:
            attract_rows.append(m)
    write_tsv_gz(
        out_dir / "attract_arabidopsis_motifs.tsv.gz",
        ["gene_id", "gene_name", "motif", "family", "len", "pubmed", "matrix_id"],
        attract_rows,
    )

    # Per-protein flat table
    flat_rows = []
    for pid, a in sorted(ann.items(), key=lambda kv: (-(kv[1]["score"] or 0), kv[0])):
        flat_rows.append(
            {
                "protein_id": pid,
                "gene_name": a["gene_name"],
                "tier": a["tier"],
                "score": f"{a['score']:.6f}" if a["score"] is not None else "",
                "ortholog_ath": a["ortholog_ath"] or "",
                "ortholog_ath_gene": a["ortholog_ath_gene"] or "",
                "source_tags": a["source_tags"],
                "rbp2go_score": (
                    f"{a['rbp2go']['rbp2go_score']:.4g}" if a.get("rbp2go") else ""
                ),
                "rbp2go_gene": a["rbp2go"]["gene_name"] if a.get("rbp2go") else "",
                "n_rgg": a["rg_motif"]["n_rgg"],
                "n_rg": a["rg_motif"]["n_rg"],
                "max_rgg_run": a["rg_motif"]["max_rgg_run"],
                "rg_motif_call": int(a["rg_motif"]["rg_motif_call"]),
                "uniprot_go_rna_binding": int(a["uniprot_go_rna_binding"]),
                "uniprot_kw_rna_binding": int(a["uniprot_kw_rna_binding"]),
                "attract_n_motifs": len(a["attract_motifs"]),
                "motif_note": a["motif_note"][:500],
                "rna_score_context": int(a["rna_score_context"]),
            }
        )
    write_tsv_gz(
        out_dir / "grape_rbp_annotations.tsv.gz",
        list(flat_rows[0].keys()) if flat_rows else ["protein_id"],
        flat_rows,
    )

    # Sidecar for UI (compact)
    sidecar = {}
    for pid, a in ann.items():
        sidecar[pid] = {
            "rbp_flag": 1,
            "sources": a["sources"],
            "motif_note": a["motif_note"][:400],
            "n_rgg": a["rg_motif"]["n_rgg"],
            "n_rg": a["rg_motif"]["n_rg"],
            "max_rgg_run": a["rg_motif"]["max_rgg_run"],
            "rg_motif_call": a["rg_motif"]["rg_motif_call"],
            "rg_motif_note": a["rg_motif"]["rg_motif_note"],
            "ortholog_ath": a["ortholog_ath"],
            "rbp2go_score": a["rbp2go"]["rbp2go_score"] if a.get("rbp2go") else None,
            "rbp2go_gene": a["rbp2go"]["gene_name"] if a.get("rbp2go") else None,
            "attract_n": len(a["attract_motifs"]),
            "attract_motifs": [
                {"motif": m["motif"], "family": m["family"]} for m in a["attract_motifs"][:8]
            ],
            "uniprot_go": a["uniprot_go_rna_binding"],
            "uniprot_kw": a["uniprot_kw_rna_binding"],
            "rna_score_context": a["rna_score_context"],
            "tier": a["tier"],
            "score": a["score"],
            "gene_name": a["gene_name"],
            "protein_name": (a["protein_name"] or "")[:120],
        }
    with gzip.open(out_dir / "annotations_by_protein.json.gz", "wt", encoding="utf-8") as fh:
        json.dump(sidecar, fh, separators=(",", ":"))

    # Browse tables (Tier A∪B priority + RG∩RBP2GO)
    def browse_row(a: dict) -> dict:
        return {
            "protein_id": a["protein_id"],
            "gene_name": a["gene_name"],
            "tier": a["tier"],
            "score": a["score"],
            "sources": a["sources"],
            "motif_note": a["motif_note"][:220],
            "n_rgg": a["rg_motif"]["n_rgg"],
            "n_rg": a["rg_motif"]["n_rg"],
            "rg_motif_call": a["rg_motif"]["rg_motif_call"],
            "rbp2go_score": a["rbp2go"]["rbp2go_score"] if a.get("rbp2go") else None,
            "ortholog_ath": a["ortholog_ath"],
        }

    tier_ab_rows = [
        browse_row(a)
        for a in sorted(ann.values(), key=lambda x: -(x["score"] or 0))
        if a["tier"] in {"A", "B"}
    ][:500]
    rg_priority = [
        browse_row(a)
        for a in sorted(ann.values(), key=lambda x: (-x["rg_motif"]["n_rgg"], -(x["score"] or 0)))
        if a["rg_motif"]["rg_motif_call"]
    ][:400]
    ui = {
        "tier_ab_rbp": tier_ab_rows,
        "rg_motif_hits": rg_priority,
        "n_annotated": n_rbp,
    }
    (out_dir / "ui_tables.json").write_text(json.dumps(ui, indent=2) + "\n", encoding="utf-8")

    summary = {
        "version": "1.9.0",
        "layer": "rbp_rna_binding_crosswalk",
        "n_atlas_proteins": n_proteins,
        "n_with_rbp_annotation": n_rbp,
        "hit_rate_proteome": round(n_rbp / n_proteins, 6) if n_proteins else 0,
        "by_source": dict(by_src),
        "by_tier": dict(tier_c),
        "n_tier_ab_with_rbp": n_tier_ab_rbp,
        "n_rg_motif_call": n_rg_call,
        "n_rbp2go_ath_entries": len(rbp2go),
        "n_rbp2go_ath_with_grape_ortholog": sum(
            1 for uid in rbp2go if ath_to_grape.get(uid)
        ),
        "n_grape_via_rbp2go": n_rbp2go_grape,
        "n_grape_via_attract": n_attract_grape,
        "n_attract_ath_uniprot_resolved": sum(1 for v in agi_map.values() if v),
        "n_attract_ath_uniprot_with_grape_ortholog": sum(
            1
            for uids in agi_map.values()
            for uid in uids
            if ath_to_grape.get(uid)
        ),
        "n_rg_and_rbp2go": n_rg_and_rbp2go,
        "rg_thresholds": {"rgg_min": rgg_min, "rg_min": rg_min},
        "species_coverage": {
            "rbp2go_species": sorted({r["species"] for r in rbp2go.values()}),
            "rbp2go_includes_vitis": False,
            "attract_n_organisms": attract_summary["n_organisms"],
            "attract_ath_genes": attract_summary["n_ath_genes"],
            "attract_ath_rows": attract_summary["n_ath_rows"],
            "attract_includes_vitis": False,
            "grape_native_uniprot_go3723": by_src.get("uniprot_go3723", 0),
            "grape_native_uniprot_kw": by_src.get("uniprot_kw_rna_binding", 0),
        },
        "attract_summary": attract_summary,
        "paths": {
            "rbp2go_tsv": "data/rbp/rbp2go_arabidopsis.tsv.gz",
            "attract_tsv": "data/rbp/attract_arabidopsis_motifs.tsv.gz",
            "grape_tsv": "data/rbp/grape_rbp_annotations.tsv.gz",
            "sidecar": "data/rbp/annotations_by_protein.json.gz",
            "ui_tables": "data/rbp/ui_tables.json",
            "method": "data/rbp/METHOD.md",
        },
        "honesty": (
            "RBP2GO 2.0 has no Vitis table — grape joins are Ath OrthoDB transfers. "
            "ATtRACT Ath motif coverage is sparse and animal-dominated. "
            "RG/RGG calls are sequence grammar thresholds, not experimental RNA binding. "
            "UniProt GO/KW are native grape annotations. Sidecar-only (no fat SQLite columns)."
        ),
    }
    (out_dir / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    method = f"""# RNA-binding annotation crosswalk (NAR layer A3)

PlantLLPS / GrapeCondensatome nucleic-acid scope layer: cross LLPS-scored grape
proteins against public RBP resources and RG/RGG sequence grammar.

## Sources

| Resource | Role | Plant / grape coverage |
|---|---|---|
| **RBP2GO 2.0** Ath CSV | Primary experimental+GO RBP catalog | Ath only (2 352 RBPs); **no Vitis** |
| **ATtRACT** | RBP-associated RNA motifs | Animal-heavy; Ath genes ≈ {attract_summary['n_ath_genes']} (sparse); **no Vitis** |
| UniProt UP000009183 | Native GO:0003723 + KW RNA-binding | Direct grape annotations |
| Grape FASTA | RG / RGG motif scan | Whole scored proteome |

Import: `scripts/import_rbp_sources.py` → `data/rbp/cache/`.  
Crosswalk: `scripts/enrich_rbp_crosswalk.py`.

## Join logic

1. Map each grape protein’s `ortholog_ath` to RBP2GO Ath `Uniprot_ID`.
2. Attach ATtRACT Ath motif rows when the Ath ortholog / gene resolves.
3. Flag native UniProt GO:0003723 and Keyword RNA-binding.
4. Scan sequences for overlapping `RGG` and `RG`; call when
   `n_RGG ≥ {rgg_min}` **or** `n_RG ≥ {rg_min}` **or** max consecutive RGG-run ≥ 2.

## Headline counts (this deposit)

- Atlas proteins: **{n_proteins}**
- Any RBP-related annotation: **{n_rbp}** ({100 * n_rbp / n_proteins:.2f}% proteome)
- Via RBP2GO Ath ortholog: **{n_rbp2go_grape}**
- Via ATtRACT Ath motifs: **{n_attract_grape}**
- RG/RGG motif call: **{n_rg_call}**
- Tier A∪B with any RBP annotation: **{n_tier_ab_rbp}**
- RG motif ∩ RBP2GO: **{n_rg_and_rbp2go}**

## Honesty

- OrthoDB groups can be broad; an Ath RBP ortholog is **not** proof the grape
  protein binds RNA.
- ATtRACT motif strings are transferred labels for triage, not grape CLIP peaks.
- RG/RGG thresholds are conventional grammar filters used to prioritize
  RNA-related LLPS candidates (stress-granule / RGG-box biology), not binders.
- No new species condensatome expansion in this layer; no Zenodo touch.
"""
    (out_dir / "METHOD.md").write_text(method, encoding="utf-8")

    readme = """# data/rbp — RNA-binding crosswalk deposits

NAR-priority nucleic-acid annotation layer (v1.9.0).

| File | Contents |
|---|---|
| `rbp2go_arabidopsis.tsv.gz` | Slim RBP2GO 2.0 Ath RBP table |
| `attract_arabidopsis_motifs.tsv.gz` | ATtRACT Ath motif rows |
| `grape_rbp_annotations.tsv.gz` | Per-grape-protein crosswalk |
| `annotations_by_protein.json.gz` | UI sidecar |
| `ui_tables.json` | Browse slices (Tier A∪B + RG hits) |
| `summary.json` | Hit rates / species coverage |
| `METHOD.md` | Methods blurb |
| `cache/` | Raw RBP_AT.csv + ATtRACT.zip (gitignored) |

Rebuild::

    python3 scripts/import_rbp_sources.py --rbp2go-csv RBP_AT.csv --fetch-attract
    python3 scripts/enrich_rbp_crosswalk.py
"""
    (out_dir / "README.md").write_text(readme, encoding="utf-8")

    # SQLite meta only
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("rbp_crosswalk", json.dumps({k: summary[k] for k in summary if k != "attract_summary"})),
    )
    con.execute(
        "INSERT INTO meta(key, value) VALUES(?, ?) "
        "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
        ("tool_version", "1.9.0"),
    )
    con.commit()
    con.close()

    print(json.dumps({k: summary[k] for k in (
        "n_atlas_proteins",
        "n_with_rbp_annotation",
        "hit_rate_proteome",
        "by_source",
        "by_tier",
        "n_grape_via_rbp2go",
        "n_grape_via_attract",
        "n_rg_motif_call",
        "n_rg_and_rbp2go",
        "species_coverage",
    )}, indent=2))
    print(f"Deposits → {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
