"""Shared HTTP helpers for enrichment scripts (polite User-Agent + retries)."""

from __future__ import annotations

import time
import urllib.error
import urllib.request
from typing import Optional


# Prefer a browser-like UA: AlphaFold file CDN returns HTTP 403 to bot-like agents.
DEFAULT_UA = (
    "Mozilla/5.0 (X11; Linux x86_64) AppleWebKit/537.36 "
    "(KHTML, like Gecko) Chrome/120.0.0.0 Safari/537.36"
)


def urlopen_bytes(
    url: str,
    *,
    timeout: float = 60.0,
    retries: int = 5,
    ua: str = DEFAULT_UA,
    accept: Optional[str] = None,
) -> bytes:
    headers = {"User-Agent": ua}
    if accept:
        headers["Accept"] = accept
    req = urllib.request.Request(url, headers=headers)
    last: Exception | None = None
    for attempt in range(retries):
        try:
            with urllib.request.urlopen(req, timeout=timeout) as resp:
                return resp.read()
        except (urllib.error.HTTPError, urllib.error.URLError, TimeoutError) as e:
            last = e
            code = getattr(e, "code", None)
            if code in (404, 400):
                raise
            time.sleep(min(32.0, 1.5 * (2**attempt)))
    assert last is not None
    raise last


def urlopen_text(url: str, **kwargs) -> str:
    return urlopen_bytes(url, **kwargs).decode("utf-8", errors="replace")
