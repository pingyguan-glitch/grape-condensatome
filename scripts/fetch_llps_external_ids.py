#!/usr/bin/env python3
"""Fetch/cache UniProt-mapped ID lists from PhaSePro, LLPSDB, DrLLPS, CD-CODE.

Writes compact redistributable ID tables under ``data/llps_crosswalk/ids/``.
Raw dumps land in ``data/llps_crosswalk/cache/`` (gitignored).

Stdlib only. Prefer dumps / UniProt-mapped tables over scraping HTML UIs.

Usage:
  python3 scripts/fetch_llps_external_ids.py
  python3 scripts/fetch_llps_external_ids.py --skip-network   # rebuild IDs from cache
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
import zipfile
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import urlopen_bytes  # noqa: E402

OUT_DIR = ROOT / "data" / "llps_crosswalk"
IDS_DIR = OUT_DIR / "ids"
CACHE_DIR = OUT_DIR / "cache"

UNIPROT_RE = re.compile(
    r"^[OPQ][0-9][A-Z0-9]{3}[0-9]$|^[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2}$"
)

SOURCES = {
    "phasepro": {
        "url": "https://phasepro.elte.hu/download_full.json",
        "cache_name": "phasepro_download_full.json",
        "citation": "Mészáros et al. 2020 Nucleic Acids Res. 48:D360–D367",
        "homepage": "https://phasepro.elte.hu/",
        "license_note": (
            "Public JSON dump of curated LLPS-driver entries (UniProt accession keys). "
            "Redistribute UniProt ID list with citation; do not claim PhaSePro curation."
        ),
    },
    "llpsdb": {
        "url": "http://bio-comp.org.cn/llpsdbv2/download/Phase_separation_unambiguous.zip",
        "cache_name": "llpsdb_Phase_separation_unambiguous.zip",
        "citation": "Li et al. 2020 Nucleic Acids Res. 48:D320–D327; LLPSDB v2 (2021)",
        "homepage": "http://bio-comp.org.cn/llpsdbv2/",
        "license_note": (
            "Free download of unambiguous phase-separation protein.xls from LLPSDB v2. "
            "We deposit UniProt accessions only (not sequences or full XLS)."
        ),
    },
    "drllps_experimental": {
        "url": "http://llps.biocuckoo.cn/download/LLPS.txt",
        "cache_name": "drllps_LLPS.txt",
        "citation": "Ning et al. 2020 Nucleic Acids Res. 48:D288–D295",
        "homepage": "http://llps.biocuckoo.cn/",
        "license_note": (
            "Download page TXT of experimentally identified LLPS-associated proteins "
            "(~9.2k rows with UniProt IDs). Free for academic use; cite DrLLPS."
        ),
    },
    "drllps_vitis": {
        "url": "http://llps.biocuckoo.cn/download/protein/Vitis_vinifera.zip",
        "cache_name": "drllps_Vitis_vinifera.zip",
        "citation": "Ning et al. 2020 Nucleic Acids Res. 48:D288–D295",
        "homepage": "http://llps.biocuckoo.cn/",
        "license_note": (
            "Species FASTA of DrLLPS ortholog-expanded Vitis entries (gene IDs, not UniProt). "
            "We deposit gene-symbol list only; raw FASTA stays in gitignored cache."
        ),
    },
    "cdcode": {
        "url": "https://edmond.mpg.de/api/access/datafile/350809",
        "cache_name": "cdcode_edmond_350809.zip",
        "citation": (
            "Rostam et al. 2023 Nat Methods (CD-CODE); Hadarovich et al. 2025 NAR "
            "(CD-CODE 2.0); Edmond dump doi:10.17617/3.WDM3PG"
        ),
        "homepage": "https://cd-code.org/release",
        "license_note": (
            "CD-CODE dumps are CC BY-NC 4.0. We deposit UniProt accession lists extracted "
            "from protein2cdcode / protein CSV for cross-annotation only — not the full dump. "
            "Commercial redistribution of CD-CODE content requires their license terms."
        ),
    },
}


def _norm_uniprot(raw: str | None) -> str | None:
    if not raw:
        return None
    tok = str(raw).strip().split()[0]
    tok = tok.split("-")[0].strip()
    if UNIPROT_RE.match(tok):
        return tok
    return None


def _write_lines(path: Path, values: set[str] | list[str], header: str) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = sorted({v for v in values if v})
    with path.open("w", encoding="utf-8") as fh:
        fh.write(f"# {header}\n")
        for v in rows:
            fh.write(v + "\n")
    return len(rows)


def _ensure_cache(name: str, url: str, *, skip_network: bool) -> Path:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    path = CACHE_DIR / name
    if path.exists() and path.stat().st_size > 0:
        return path
    if skip_network:
        raise FileNotFoundError(f"Missing cache {path}; run without --skip-network")
    print(f"GET {url} → {path.name}", flush=True)
    data = urlopen_bytes(url, timeout=180.0, retries=5)
    path.write_bytes(data)
    return path


def parse_phasepro(path: Path) -> set[str]:
    data = json.loads(path.read_text(encoding="utf-8"))
    ids: set[str] = set()
    if isinstance(data, dict):
        for k, v in data.items():
            u = _norm_uniprot(k)
            if u:
                ids.add(u)
            if isinstance(v, dict):
                u2 = _norm_uniprot(v.get("accession") or v.get("uniprot") or "")
                if u2:
                    ids.add(u2)
    return ids


def _xlsx_shared_strings(zf: zipfile.ZipFile) -> list[str]:
    raw = zf.read("xl/sharedStrings.xml")
    root = ET.fromstring(raw)
    ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
    out: list[str] = []
    for si in root.findall("m:si", ns):
        texts = [t.text or "" for t in si.findall(".//m:t", ns)]
        out.append("".join(texts))
    return out


def _xlsx_cell_value(cell: ET.Element, shared: list[str], ns: dict) -> str | None:
    t = cell.attrib.get("t")
    v = cell.find("m:v", ns)
    if v is None or v.text is None:
        return None
    if t == "s":
        return shared[int(v.text)]
    return v.text


def parse_llpsdb_zip(path: Path) -> set[str]:
    """Extract UniProt IDs from protein.xls (xlsx bytes) inside the LLPSDB zip."""
    ids: set[str] = set()
    with zipfile.ZipFile(path) as outer:
        xls_name = None
        for name in outer.namelist():
            if name.lower().endswith("protein.xls") or name.lower().endswith("protein.xlsx"):
                xls_name = name
                break
        if not xls_name:
            raise RuntimeError(f"No protein.xls in {path}")
        xls_bytes = outer.read(xls_name)

    with zipfile.ZipFile(io.BytesIO(xls_bytes)) as zf:
        shared = _xlsx_shared_strings(zf)
        ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        sheet = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        rows = sheet.findall("m:sheetData/m:row", ns)
        if not rows:
            return ids
        # Header → column letter for UniProt ID
        header_cells = rows[0].findall("m:c", ns)
        col_uniprot = None
        for c in header_cells:
            ref = c.attrib.get("r", "")
            col = "".join(ch for ch in ref if ch.isalpha())
            val = _xlsx_cell_value(c, shared, ns) or ""
            if "uniprot" in val.lower():
                col_uniprot = col
                break
        if not col_uniprot:
            col_uniprot = "F"  # LLPSDB v2 layout
        for row in rows[1:]:
            vals: dict[str, str | None] = {}
            for c in row.findall("m:c", ns):
                ref = c.attrib.get("r", "")
                col = "".join(ch for ch in ref if ch.isalpha())
                vals[col] = _xlsx_cell_value(c, shared, ns)
            raw = vals.get(col_uniprot)
            if not raw:
                continue
            for tok in re.split(r"[\s;,/|]+", str(raw).strip()):
                u = _norm_uniprot(tok)
                if u:
                    ids.add(u)
    return ids


def parse_drllps_experimental(path: Path) -> set[str]:
    ids: set[str] = set()
    with path.open(encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # Header: DrLLPS ID / UniProt ID / ...
        for row in reader:
            raw = row.get("UniProt ID") or row.get("UniProt") or ""
            if not raw:
                # positional fallback
                vals = list(row.values())
                raw = vals[1] if len(vals) > 1 else ""
            u = _norm_uniprot(raw)
            if u:
                ids.add(u)
    return ids


def parse_drllps_vitis_genes(path: Path) -> set[str]:
    genes: set[str] = set()
    with zipfile.ZipFile(path) as zf:
        names = [n for n in zf.namelist() if not n.endswith("/")]
        if not names:
            raise RuntimeError(f"Empty DrLLPS Vitis zip: {path}")
        text = zf.read(names[0]).decode("utf-8", errors="replace")
    for line in text.splitlines():
        if not line.startswith(">"):
            continue
        # >LLPS-Viv-0000|VIT_04s0023g03300|Vitis vinifera
        parts = line[1:].strip().split("|")
        if len(parts) >= 2 and parts[1].strip():
            genes.add(parts[1].strip())
    return genes


def parse_cdcode_zip(path: Path) -> set[str]:
    ids: set[str] = set()
    with zipfile.ZipFile(path) as zf:
        # Prefer protein2cdcode_*.tsv
        tsv_names = [n for n in zf.namelist() if "protein2cdcode" in n.lower() and n.endswith(".tsv")]
        csv_names = [
            n
            for n in zf.namelist()
            if n.lower().startswith("protein_") and n.endswith(".csv") and not n.startswith("__")
        ]
        if tsv_names:
            text = zf.read(tsv_names[0]).decode("utf-8", errors="replace")
            reader = csv.DictReader(io.StringIO(text), delimiter="\t")
            for row in reader:
                u = _norm_uniprot(row.get("uniprotkb_ac") or row.get("uniprot_id") or "")
                if u:
                    ids.add(u)
        if csv_names:
            text = zf.read(csv_names[0]).decode("utf-8", errors="replace")
            reader = csv.DictReader(io.StringIO(text))
            for row in reader:
                u = _norm_uniprot(row.get("uniprot_id") or row.get("uniprotkb_id") or "")
                if u:
                    ids.add(u)
                # uniprotkb_id sometimes like Q9U302_CAEEL
                raw = row.get("uniprotkb_id") or ""
                if "_" in raw:
                    u2 = _norm_uniprot(raw.split("_")[0])
                    if u2:
                        ids.add(u2)
    return ids


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--skip-network",
        action="store_true",
        help="Rebuild ID tables from existing cache/ only",
    )
    args = ap.parse_args()

    IDS_DIR.mkdir(parents=True, exist_ok=True)
    CACHE_DIR.mkdir(parents=True, exist_ok=True)

    provenance: dict = {
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "sources": {},
        "counts": {},
    }

    # PhaSePro
    pp_path = _ensure_cache(
        SOURCES["phasepro"]["cache_name"],
        SOURCES["phasepro"]["url"],
        skip_network=args.skip_network,
    )
    phasepro = parse_phasepro(pp_path)
    n = _write_lines(
        IDS_DIR / "phasepro_uniprot.txt",
        phasepro,
        "PhaSePro UniProt accessions (from download_full.json keys/accession)",
    )
    provenance["sources"]["phasepro"] = {**SOURCES["phasepro"], "cache": str(pp_path.relative_to(ROOT))}
    provenance["counts"]["phasepro_uniprot"] = n

    # LLPSDB
    ll_path = _ensure_cache(
        SOURCES["llpsdb"]["cache_name"],
        SOURCES["llpsdb"]["url"],
        skip_network=args.skip_network,
    )
    llpsdb = parse_llpsdb_zip(ll_path)
    n = _write_lines(
        IDS_DIR / "llpsdb_uniprot.txt",
        llpsdb,
        "LLPSDB v2 unambiguous phase-separation UniProt IDs (protein.xls)",
    )
    provenance["sources"]["llpsdb"] = {**SOURCES["llpsdb"], "cache": str(ll_path.relative_to(ROOT))}
    provenance["counts"]["llpsdb_uniprot"] = n

    # DrLLPS experimental
    dr_path = _ensure_cache(
        SOURCES["drllps_experimental"]["cache_name"],
        SOURCES["drllps_experimental"]["url"],
        skip_network=args.skip_network,
    )
    drllps = parse_drllps_experimental(dr_path)
    n = _write_lines(
        IDS_DIR / "drllps_experimental_uniprot.txt",
        drllps,
        "DrLLPS experimentally identified LLPS-associated UniProt IDs (LLPS.txt)",
    )
    provenance["sources"]["drllps_experimental"] = {
        **SOURCES["drllps_experimental"],
        "cache": str(dr_path.relative_to(ROOT)),
    }
    provenance["counts"]["drllps_experimental_uniprot"] = n

    # DrLLPS Vitis gene symbols
    viv_path = _ensure_cache(
        SOURCES["drllps_vitis"]["cache_name"],
        SOURCES["drllps_vitis"]["url"],
        skip_network=args.skip_network,
    )
    viv_genes = parse_drllps_vitis_genes(viv_path)
    n = _write_lines(
        IDS_DIR / "drllps_vitis_genes.txt",
        viv_genes,
        "DrLLPS Vitis vinifera ortholog-expansion gene IDs (VIT_…)",
    )
    provenance["sources"]["drllps_vitis"] = {
        **SOURCES["drllps_vitis"],
        "cache": str(viv_path.relative_to(ROOT)),
    }
    provenance["counts"]["drllps_vitis_genes"] = n

    # CD-CODE
    cd_path = _ensure_cache(
        SOURCES["cdcode"]["cache_name"],
        SOURCES["cdcode"]["url"],
        skip_network=args.skip_network,
    )
    cdcode = parse_cdcode_zip(cd_path)
    n = _write_lines(
        IDS_DIR / "cdcode_uniprot.txt",
        cdcode,
        "CD-CODE UniProt accessions (Edmond protein2cdcode / protein CSV; CC BY-NC)",
    )
    provenance["sources"]["cdcode"] = {**SOURCES["cdcode"], "cache": str(cd_path.relative_to(ROOT))}
    provenance["counts"]["cdcode_uniprot"] = n

    # Union / intersections for quick sanity
    provenance["counts"]["union_uniprot_phasepro_llpsdb_drllps_cdcode"] = len(
        phasepro | llpsdb | drllps | cdcode
    )
    meta_path = IDS_DIR / "sources.json"
    meta_path.write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(provenance["counts"], indent=2))
    print(f"Wrote {meta_path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
