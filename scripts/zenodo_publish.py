#!/usr/bin/env python3
"""Publish a prepared PlantLLPS / GrapeCondensatome archive to Zenodo via REST API.

Requires a personal access token (never commit tokens):

    export ZENODO_ACCESS_TOKEN=…   # or ZENODO_TOKEN
    python3 scripts/zenodo_publish.py --archive dist/grape-condensatome-zenodo-v1.11.0.zip

Dry-run (no network writes beyond auth check):

    python3 scripts/zenodo_publish.py --archive … --dry-run

After publish, paste the printed DOI into CITATION.cff (see docs/ZENODO.md).
This script never invents a DOI.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DEFAULT_API = "https://zenodo.org/api"


def _token() -> str | None:
    for key in ("ZENODO_ACCESS_TOKEN", "ZENODO_TOKEN"):
        val = os.environ.get(key, "").strip()
        if val:
            return val
    return None


def _read_version() -> str:
    path = ROOT / "VERSION"
    if path.is_file():
        return path.read_text(encoding="utf-8").strip().lstrip("v") or "1.11.0"
    return "1.11.0"


def _api_request(
    method: str,
    url: str,
    token: str,
    *,
    data: bytes | None = None,
    content_type: str | None = None,
    timeout: float = 120.0,
) -> dict:
    headers = {"Authorization": f"Bearer {token}", "Accept": "application/json"}
    if content_type:
        headers["Content-Type"] = content_type
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            body = resp.read().decode("utf-8")
            return json.loads(body) if body else {}
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"Zenodo API {method} {url} → HTTP {exc.code}: {detail}") from exc


def _metadata(version: str) -> dict:
    abstract = (
        "PlantLLPS is a plant-aware sequence scorer for liquid–liquid phase separation "
        "propensity. GrapeCondensatome is the PlantLLPS-scored UniProt UP000009183 "
        "Vitis vinifera proteome atlas with confidence tiers and UniProt annotation joins. "
        "This deposit packages PlantLLPS v13 weights, the grape atlas SQLite/TSV exports, "
        "a paste-FASTA scoring-server snapshot, and the sample FASTA used by /score "
        f"(software tip v{version}). Gold-panel AUCs are partially in-sample; production "
        "disorder is a composition proxy. Do not invent alternate DOIs — cite this record."
    )
    return {
        "metadata": {
            "title": f"PlantLLPS and GrapeCondensatome (v{version})",
            "upload_type": "software",
            "description": abstract,
            "creators": [
                {
                    "name": "Guan, Pingyin",
                    "affiliation": (
                        "College of Horticulture, China Agricultural University, "
                        "Beijing 100193, China"
                    ),
                    "orcid": "0000-0003-4890-3412",
                }
            ],
            "version": version,
            "license": "mit",
            "keywords": [
                "liquid-liquid phase separation",
                "condensatome",
                "Vitis vinifera",
                "PlantLLPS",
                "web server",
            ],
            "related_identifiers": [
                {
                    "identifier": "https://github.com/pingyguan-glitch/grape-condensatome",
                    "relation": "isSupplementTo",
                    "resource_type": "software",
                    "scheme": "url",
                },
                {
                    "identifier": (
                        "https://github.com/pingyguan-glitch/grape-condensatome"
                        f"/releases/tag/v{version}"
                    ),
                    "relation": "isIdenticalTo",
                    "resource_type": "software",
                    "scheme": "url",
                },
            ],
            "access_right": "open",
        }
    }


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--archive",
        type=Path,
        required=True,
        help="Path to deposit .zip or .tar.gz (from export_zenodo_bundle.py)",
    )
    ap.add_argument("--version", default=None, help="Version label (default: VERSION file)")
    ap.add_argument(
        "--api",
        default=os.environ.get("ZENODO_API", DEFAULT_API),
        help="Zenodo API root (default: https://zenodo.org/api; sandbox: https://sandbox.zenodo.org/api)",
    )
    ap.add_argument(
        "--publish",
        action="store_true",
        help="Publish the deposition after upload (mints the DOI). Without this flag, only a draft is created.",
    )
    ap.add_argument(
        "--dry-run",
        action="store_true",
        help="Validate token + archive; print planned metadata; do not create a deposition",
    )
    args = ap.parse_args(argv)

    archive = args.archive.expanduser().resolve()
    if not archive.is_file():
        raise SystemExit(f"Missing archive: {archive}")

    version = (args.version or _read_version()).lstrip("v")
    token = _token()
    if not token:
        print(
            "No Zenodo token in the environment.\n"
            "  1. Create a personal access token at https://zenodo.org/account/settings/applications/\n"
            "     (scopes: deposit:actions, deposit:write)\n"
            "  2. export ZENODO_ACCESS_TOKEN=…   # or ZENODO_TOKEN\n"
            "  3. Re-run this script, or upload the archive manually (docs/ZENODO.md).\n"
            f"Archive ready: {archive} ({archive.stat().st_size / (1 << 20):.1f} MiB)",
            file=sys.stderr,
        )
        return 2

    meta = _metadata(version)
    print(f"Archive: {archive} ({archive.stat().st_size / (1 << 20):.1f} MiB)")
    print(f"API:     {args.api}")
    print(f"Version: {version}")
    print(f"Publish: {args.publish}")
    print("Metadata preview:")
    print(json.dumps(meta, indent=2))

    if args.dry_run:
        # Auth check only
        _api_request("GET", f"{args.api.rstrip('/')}/deposit/depositions?size=1", token)
        print("Dry-run OK (token accepted). No deposition created.")
        return 0

    base = args.api.rstrip("/")
    dep = _api_request(
        "POST",
        f"{base}/deposit/depositions",
        token,
        data=b"{}",
        content_type="application/json",
    )
    dep_id = dep["id"]
    bucket = dep["links"]["bucket"]
    print(f"Created draft deposition id={dep_id}")

    # Upload file into the bucket
    upload_url = f"{bucket}/{archive.name}"
    with archive.open("rb") as fh:
        raw = fh.read()
    req = urllib.request.Request(
        upload_url,
        data=raw,
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/octet-stream",
            "Accept": "application/json",
        },
        method="PUT",
    )
    try:
        with urllib.request.urlopen(req, timeout=600.0) as resp:
            _ = resp.read()
    except urllib.error.HTTPError as exc:
        detail = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"Upload failed HTTP {exc.code}: {detail}") from exc
    print(f"Uploaded {archive.name}")

    updated = _api_request(
        "PUT",
        f"{base}/deposit/depositions/{dep_id}",
        token,
        data=json.dumps(meta).encode("utf-8"),
        content_type="application/json",
    )
    print(f"Metadata set. Draft URL: {updated['links'].get('html', updated['links'].get('latest_html'))}")

    if not args.publish:
        print(
            "Draft only (pass --publish to mint the DOI).\n"
            "Open the draft URL in the browser, review, then Publish — or re-run with --publish."
        )
        return 0

    published = _api_request(
        "POST",
        f"{base}/deposit/depositions/{dep_id}/actions/publish",
        token,
    )
    doi = published.get("doi") or published.get("metadata", {}).get("prereserve_doi", {}).get("doi")
    concept = published.get("conceptdoi")
    html = published["links"].get("html") or published["links"].get("latest_html")
    print("Published.")
    print(f"Version DOI:  {doi}")
    if concept:
        print(f"Concept DOI:  {concept}")
    print(f"Record URL:   {html}")
    print(
        "Next: paste the Version DOI into CITATION.cff identifiers "
        "(see docs/ZENODO.md). Do not invent a different DOI."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
