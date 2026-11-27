#!/usr/bin/env python3
"""Import / cache public RBP resources used by the NAR nucleic-acid crosswalk.

Sources
-------
* **RBP2GO 2.0** Arabidopsis thaliana RBP dataset (CSV via the Shiny
  “References & Data → Download RBP Dataset” control). *Vitis* is **not** a
  RBP2GO species — plant coverage is Ath-only; grape joins use OrthoDB.
* **ATtRACT** RBP + motif table (ZIP from https://attract.cnic.es/).
  Organism mix is animal-heavy; Ath motif rows are sparse (documented gap).

Caches land under ``data/rbp/cache/`` (gitignored bulk). Slim deposited tables
are written by ``scripts/enrich_rbp_crosswalk.py``.

Examples::

    # Prefer a previously downloaded Ath CSV (Playwright / browser export)
    python3 scripts/import_rbp_sources.py \\
        --rbp2go-csv /path/to/RBP_AT.csv

    # Re-fetch ATtRACT; attempt Playwright download of RBP2GO when CSV missing
    python3 scripts/import_rbp_sources.py --fetch-attract --fetch-rbp2go
"""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from _http_util import urlopen_bytes  # noqa: E402

DEFAULT_DIR = ROOT / "data" / "rbp"
CACHE_DIR = DEFAULT_DIR / "cache"
ATTRACT_URL = "https://attract.cnic.es/attract/static/ATtRACT.zip"
RBP2GO_HOME = "https://rbp2gov2.dkfz.de/"


def _skip_hash_comments(fh) -> None:
    while True:
        pos = fh.tell()
        line = fh.readline()
        if not line:
            return
        if not line.startswith("#"):
            fh.seek(pos)
            return


def validate_rbp2go_csv(path: Path) -> dict:
    with path.open(newline="", encoding="utf-8", errors="replace") as fh:
        _skip_hash_comments(fh)
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SystemExit(f"empty RBP2GO CSV: {path}")
    need = {"Uniprot_ID", "Gene_Name", "RBP2GO_Score", "Species"}
    missing = need - set(rows[0].keys())
    if missing:
        raise SystemExit(f"RBP2GO CSV missing columns {sorted(missing)}")
    species = {r.get("Species", "") for r in rows}
    return {
        "n_rows": len(rows),
        "species": sorted(s for s in species if s),
        "path": str(path),
    }


def fetch_attract(cache_dir: Path, force: bool = False) -> Path:
    dest = cache_dir / "ATtRACT.zip"
    if dest.exists() and dest.stat().st_size > 1000 and not force:
        print(f"Using ATtRACT cache → {dest}", flush=True)
        return dest
    cache_dir.mkdir(parents=True, exist_ok=True)
    print(f"Downloading ATtRACT → {dest}", flush=True)
    dest.write_bytes(urlopen_bytes(ATTRACT_URL, timeout=120))
    return dest


def fetch_rbp2go_playwright(cache_dir: Path) -> Path:
    """Best-effort headless download of Ath RBP_AT.csv from RBP2GO 2.0."""
    try:
        from playwright.sync_api import sync_playwright
    except ImportError as e:
        raise SystemExit(
            "Playwright not installed. Pass --rbp2go-csv PATH to a browser-exported "
            f"RBP_AT.csv from {RBP2GO_HOME} (References & Data → Arabidopsis thaliana), "
            f"or `pip install playwright && playwright install chromium`. ({e})"
        ) from e

    cache_dir.mkdir(parents=True, exist_ok=True)
    out = cache_dir / "RBP_AT.csv"
    print(f"Playwright download from {RBP2GO_HOME} …", flush=True)
    with sync_playwright() as p:
        browser = p.chromium.launch(headless=True)
        context = browser.new_context(accept_downloads=True)
        page = context.new_page()
        page.goto(RBP2GO_HOME, wait_until="domcontentloaded", timeout=120_000)
        # Open References & Data tab then click Ath download.
        # Selectors are resilient to Shiny re-renders.
        try:
            page.get_by_text("References & Data", exact=False).first.click(timeout=30_000)
        except Exception:
            page.get_by_role("link", name="References & Data").click(timeout=30_000)
        page.wait_for_timeout(1500)
        with page.expect_download(timeout=120_000) as dl_info:
            # Prefer an explicit Arabidopsis download control
            loc = page.get_by_text("Download RBP Dataset", exact=False)
            # Click the Ath-associated control: walk nearby text
            ath = page.get_by_text("Arabidopsis thaliana", exact=False)
            if ath.count():
                # click download near Ath section — fall back to first matching button
                try:
                    ath.first.locator("xpath=ancestor::*[.//a or .//button][1]").get_by_text(
                        "Download RBP Dataset", exact=False
                    ).first.click(timeout=15_000)
                except Exception:
                    loc.first.click(timeout=30_000)
            else:
                loc.first.click(timeout=30_000)
        download = dl_info.value
        download.save_as(str(out))
        browser.close()
    info = validate_rbp2go_csv(out)
    print(f"RBP2GO Ath CSV: {info['n_rows']} rows → {out}", flush=True)
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out-dir", type=Path, default=DEFAULT_DIR)
    ap.add_argument("--rbp2go-csv", type=Path, default=None, help="Pre-downloaded RBP_AT.csv")
    ap.add_argument("--fetch-attract", action="store_true", help="Download ATtRACT.zip")
    ap.add_argument(
        "--fetch-rbp2go",
        action="store_true",
        help="Attempt Playwright download of Ath RBP2GO CSV",
    )
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    cache = args.out_dir / "cache"
    cache.mkdir(parents=True, exist_ok=True)

    attract = None
    if args.fetch_attract or not (cache / "ATtRACT.zip").exists():
        attract = fetch_attract(cache, force=args.force)
    else:
        attract = cache / "ATtRACT.zip"
        print(f"Using ATtRACT cache → {attract}", flush=True)

    rbp2go = None
    if args.rbp2go_csv:
        src = args.rbp2go_csv.expanduser().resolve()
        if not src.exists():
            raise SystemExit(f"--rbp2go-csv not found: {src}")
        dest = cache / "RBP_AT.csv"
        if src.resolve() != dest.resolve():
            shutil.copy2(src, dest)
        rbp2go = dest
        print(validate_rbp2go_csv(rbp2go), flush=True)
    elif (cache / "RBP_AT.csv").exists() and not args.force:
        rbp2go = cache / "RBP_AT.csv"
        print(f"Using RBP2GO cache → {rbp2go}", flush=True)
        print(validate_rbp2go_csv(rbp2go), flush=True)
    elif args.fetch_rbp2go:
        rbp2go = fetch_rbp2go_playwright(cache)
    else:
        print(
            "WARNING: no RBP2GO Ath CSV in cache. Re-run with --rbp2go-csv or --fetch-rbp2go.",
            flush=True,
        )

    manifest = {
        "attract_zip": str(attract) if attract else None,
        "rbp2go_csv": str(rbp2go) if rbp2go else None,
        "rbp2go_home": RBP2GO_HOME,
        "attract_url": ATTRACT_URL,
        "note": (
            "RBP2GO has no Vitis species table; Ath CSV is the plant-aware source. "
            "ATtRACT Ath motif coverage is sparse — see enrich_rbp_crosswalk.py honesty."
        ),
    }
    (cache / "import_manifest.json").write_text(
        __import__("json").dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print("Wrote", cache / "import_manifest.json", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
