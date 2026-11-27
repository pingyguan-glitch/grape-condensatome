"""Stdlib HTTP server for browsing GrapeCondensatome SQLite + on-the-fly scoring."""

from __future__ import annotations

import html
import json
import re
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple
from urllib.parse import parse_qs, quote, unquote, urlparse

from .. import __version__
from ..db import (
    connect,
    get_meta,
    get_protein,
    list_by_tier,
    search_proteins,
    summarize,
    top_proteins,
)
from ..model import ensure_trained_model

WEB_ROOT = Path(__file__).resolve().parent
STATIC_ROOT = WEB_ROOT / "static"
PACKAGE_ROOT = WEB_ROOT.parent.parent
DEFAULT_DB = PACKAGE_ROOT / "data" / "grape_condensatome.sqlite"
DEFAULT_HOST = "0.0.0.0"
DEFAULT_PORT = 8765
PAGE_SIZE = 40


def _esc(s: Any) -> str:
    return html.escape("" if s is None else str(s))


def _fmt_score(x: Any) -> str:
    try:
        return f"{float(x):.4f}"
    except (TypeError, ValueError):
        return "—"


def _fmt_pct(x: Any) -> str:
    try:
        return f"{float(x):.2f}"
    except (TypeError, ValueError):
        return "—"


def _layout(title: str, body: str, *, query: str = "") -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>{_esc(title)} · GrapeCondensatome</title>
<link rel="stylesheet" href="/static/style.css"/>
</head>
<body>
<div class="wrap">
  <header class="site-header">
    <div class="brand">
      <h1>GrapeCondensatome</h1>
      <p>Browse PlantLLPS scores for the grapevine reference proteome (UniProt UP000009183).</p>
    </div>
    <nav class="nav">
      <a href="/">Home</a>
      <a href="/tiers/A">Tier A</a>
      <a href="/tiers/B">Tier B</a>
      <a href="/tiers/C">Tier C</a>
      <a href="/tiers/D">Tier D</a>
      <a href="/score">Score FASTA</a>
    </nav>
    <!-- Full /compare + /hooks live on the Node UI (server.js). -->
  </header>
  <div class="banner">
    PlantLLPS hybrid model (v13; hard disordered negatives). Disorder defaults to a composition proxy;
    metapredict + AlphaFold pLDDT where joined; Ath/rice OrthoDB orthologs; gold PSAP/LLPhyScore + proteome PLAAC-approx.
    Treat ranks as hypotheses. Tool v{_esc(__version__)}.
  </div>
  <form class="inline card" method="get" action="/search" style="margin-bottom:1.25rem">
    <label for="q" class="muted">Search UniProt ID, gene name, GO, or keyword</label>
    <input id="q" type="text" name="q" value="{_esc(query)}" placeholder="e.g. F6HET4, MAP65, VIT_19s0090g01230" autofocus/>
    <button type="submit">Search</button>
  </form>
  {body}
  <footer class="site-footer">
    PlantLLPS + GrapeCondensatome prototype · data from UniProt UP000009183 when scored ·
    run with <span class="mono">python3 -m plantllps serve</span>
  </footer>
</div>
</body>
</html>"""


def _protein_table(rows: List[Dict[str, Any]], *, empty: str = "No proteins found.") -> str:
    if not rows:
        return f'<p class="muted">{_esc(empty)}</p>'
    cells = []
    for r in rows:
        pid = _esc(r["protein_id"])
        tier = _esc(r.get("tier", "?"))
        gene = _esc(r.get("gene_name") or "—")
        pname = _esc((r.get("protein_name") or "")[:72] or "—")
        tpm_v = r.get("expr_mean_tpm")
        try:
            tpm = f"{float(tpm_v):.2f}" if tpm_v is not None else "—"
        except (TypeError, ValueError):
            tpm = "—"
        ath = _esc(r.get("ortholog_ath_gene") or r.get("ortholog_ath") or "—")
        cells.append(
            f"<tr>"
            f'<td class="mono"><a href="/protein/{pid}">{pid}</a></td>'
            f'<td class="mono">{gene}</td>'
            f'<td><span class="tier {tier}">{tier}</span></td>'
            f"<td>{_fmt_score(r.get('score'))}</td>"
            f"<td>{_fmt_pct(r.get('percentile'))}</td>"
            f"<td>{_fmt_score(r.get('map65_tail_like'))}</td>"
            f"<td>{_fmt_score(r.get('frac_disorder_proxy'))}</td>"
            f"<td>{_fmt_score(r.get('frac_disorder_metapredict'))}</td>"
            f"<td>{_fmt_score(r.get('af_frac_very_low'))}</td>"
            f"<td>{_fmt_score(r.get('plaac_nllr_approx'))}</td>"
            f'<td class="mono">{ath}</td>'
            f"<td>{_esc(tpm)}</td>"
            f"<td class='muted'>{pname}</td>"
            f"</tr>"
        )
    return f"""
    <div class="card" style="overflow-x:auto">
      <table>
        <thead>
          <tr>
            <th>Protein</th><th>Gene</th><th>Tier</th><th>PlantLLPS</th><th>Percentile</th>
            <th>MAP65-like</th><th>Disorder proxy</th><th>metapredict</th><th>AF pLDDT&lt;50</th>
            <th>PLAAC≈</th><th>Ath ortholog</th><th>Mean TPM</th><th>Protein name</th>
          </tr>
        </thead>
        <tbody>
          {''.join(cells)}
        </tbody>
      </table>
    </div>"""


def _pager(base: str, page: int, total: int, page_size: int = PAGE_SIZE) -> str:
    n_pages = max(1, (total + page_size - 1) // page_size)
    prev_disabled = page <= 1
    next_disabled = page >= n_pages
    sep = "&" if "?" in base else "?"
    prev = f'{base}{sep}page={page-1}' if not prev_disabled else "#"
    nxt = f'{base}{sep}page={page+1}' if not next_disabled else "#"
    prev_attr = ' aria-disabled="true" tabindex="-1"' if prev_disabled else ""
    next_attr = ' aria-disabled="true" tabindex="-1"' if next_disabled else ""
    return f"""
    <div class="pager">
      <a class="btn secondary" href="{_esc(prev)}"{prev_attr}>← Prev</a>
      <span class="muted">Page {page} / {n_pages} · {total} proteins</span>
      <a class="btn secondary" href="{_esc(nxt)}"{next_attr}>Next →</a>
    </div>"""


class AppState:
    def __init__(self, db_path: Path, model_path: Optional[Path] = None):
        self.db_path = Path(db_path)
        self.model_path = model_path
        self._model = None

    def conn(self):
        return connect(self.db_path)

    def model(self):
        if self._model is None:
            if self.model_path and Path(self.model_path).exists():
                from ..model import PlantLLPSModel

                self._model = PlantLLPSModel.load(self.model_path)
            else:
                self._model = ensure_trained_model(PACKAGE_ROOT)
        return self._model


def make_handler(state: AppState) -> Callable:
    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt: str, *args) -> None:  # quieter default
            sys_stderr_write = __import__("sys").stderr.write
            sys_stderr_write("%s - %s\n" % (self.address_string(), fmt % args))

        def _send(self, code: int, body: bytes, content_type: str) -> None:
            self.send_response(code)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(body)))
            self.send_header("Cache-Control", "no-store")
            self.end_headers()
            self.wfile.write(body)

        def _html(self, code: int, title: str, body: str, query: str = "") -> None:
            page = _layout(title, body, query=query)
            self._send(code, page.encode("utf-8"), "text/html; charset=utf-8")

        def _redirect(self, loc: str) -> None:
            self.send_response(302)
            self.send_header("Location", loc)
            self.end_headers()

        def do_GET(self) -> None:
            try:
                self._route_get()
            except Exception:
                tb = traceback.format_exc()
                self._html(500, "Error", f'<div class="card"><pre class="err">{_esc(tb)}</pre></div>')

        def do_POST(self) -> None:
            try:
                self._route_post()
            except Exception:
                tb = traceback.format_exc()
                self._html(500, "Error", f'<div class="card"><pre class="err">{_esc(tb)}</pre></div>')

        def _route_get(self) -> None:
            parsed = urlparse(self.path)
            path = unquote(parsed.path)
            qs = parse_qs(parsed.query)

            if path.startswith("/static/"):
                return self._static(path[len("/static/") :])

            if path in {"/", "/index.html"}:
                return self._home()
            if path == "/search":
                return self._search(qs)
            m = re.match(r"^/tiers/([A-Da-d])/?$", path)
            if m:
                return self._tier(m.group(1).upper(), qs)
            m = re.match(r"^/protein/([^/]+)/?$", path)
            if m:
                return self._protein(m.group(1))
            if path == "/score":
                return self._score_form()
            if path == "/api/health":
                body = json.dumps({"ok": True, "db": str(state.db_path), "version": __version__}).encode()
                return self._send(200, body, "application/json")
            if path == "/api/search":
                return self._api_search(qs)

            self._html(404, "Not found", '<div class="card"><p>Page not found.</p><p><a href="/">Home</a></p></div>')

        def _route_post(self) -> None:
            parsed = urlparse(self.path)
            path = unquote(parsed.path)
            length = int(self.headers.get("Content-Length", 0))
            raw = self.rfile.read(length).decode("utf-8", errors="replace")
            form = parse_qs(raw)
            if path == "/score":
                return self._score_submit(form)
            self._html(404, "Not found", '<div class="card"><p>Unknown POST.</p></div>')

        def _static(self, rel: str) -> None:
            rel = rel.lstrip("/").replace("..", "")
            fp = (STATIC_ROOT / rel).resolve()
            if not str(fp).startswith(str(STATIC_ROOT.resolve())) or not fp.is_file():
                self._send(404, b"not found", "text/plain")
                return
            data = fp.read_bytes()
            ctype = "text/css" if fp.suffix == ".css" else "application/octet-stream"
            self._send(200, data, ctype)

        def _home(self) -> None:
            conn = state.conn()
            try:
                summary = summarize(conn)
                meta = get_meta(conn)
                top = top_proteins(conn, 15)
            finally:
                conn.close()
            tiers = summary.get("tiers") or {}
            stats = "".join(
                f'<div class="card stat"><div class="n">{_esc(tiers.get(t, 0))}</div>'
                f'<div class="l"><a href="/tiers/{t}">Tier {t}</a></div></div>'
                for t in ("A", "B", "C", "D")
            )
            body = f"""
            <div class="grid stats">
              <div class="card stat"><div class="n">{summary['n_proteins']}</div><div class="l">Proteins in DB</div></div>
              {stats}
            </div>
            <div class="card" style="margin-bottom:1rem">
              <h2>Database</h2>
              <p class="muted mono">{_esc(state.db_path)}</p>
              <p>Status: <strong>{_esc(meta.get('status', 'unknown'))}</strong>
                 · Proteome: {_esc(meta.get('proteome', '—'))}
                 · n_scored meta: {_esc(meta.get('n_scored', '—'))}</p>
              <p class="muted">{_esc(meta.get('disclaimer', ''))}</p>
            </div>
            <h2 style="font-family:var(--font-display)">Top scores</h2>
            {_protein_table(top)}
            """
            self._html(200, "Home", body)

        def _search(self, qs: Dict[str, List[str]]) -> None:
            q = (qs.get("q") or [""])[0]
            page = max(1, int((qs.get("page") or ["1"])[0]))
            offset = (page - 1) * PAGE_SIZE
            conn = state.conn()
            try:
                rows, total = search_proteins(conn, q, limit=PAGE_SIZE, offset=offset)
            finally:
                conn.close()
            base = f"/search?q={quote(q)}"
            body = f"""
            <div class="card"><h2>Search results</h2>
            <p class="muted">Query: <span class="mono">{_esc(q)}</span> · {total} hit(s)</p></div>
            {_protein_table(rows)}
            {_pager(base, page, total) if total > PAGE_SIZE else ''}
            """
            self._html(200, f"Search: {q}", body, query=q)

        def _api_search(self, qs: Dict[str, List[str]]) -> None:
            q = (qs.get("q") or [""])[0]
            conn = state.conn()
            try:
                rows, total = search_proteins(conn, q, limit=20, offset=0)
            finally:
                conn.close()
            payload = {
                "query": q,
                "total": total,
                "results": [
                    {
                        "protein_id": r["protein_id"],
                        "score": r.get("score"),
                        "tier": r.get("tier"),
                        "percentile": r.get("percentile"),
                        "map65_tail_like": r.get("map65_tail_like"),
                    }
                    for r in rows
                ],
            }
            self._send(200, json.dumps(payload).encode(), "application/json")

        def _tier(self, tier: str, qs: Dict[str, List[str]]) -> None:
            page = max(1, int((qs.get("page") or ["1"])[0]))
            offset = (page - 1) * PAGE_SIZE
            conn = state.conn()
            try:
                rows, total = list_by_tier(conn, tier, limit=PAGE_SIZE, offset=offset)
            finally:
                conn.close()
            body = f"""
            <div class="card">
              <h2>Tier <span class="tier {tier}">{tier}</span></h2>
              <p class="muted">{total} proteins · sorted by PlantLLPS score</p>
            </div>
            {_protein_table(rows)}
            {_pager(f'/tiers/{tier}', page, total)}
            """
            self._html(200, f"Tier {tier}", body)

        def _protein(self, protein_id: str) -> None:
            conn = state.conn()
            try:
                p = get_protein(conn, protein_id)
            finally:
                conn.close()
            if not p:
                self._html(
                    404,
                    "Missing",
                    f'<div class="card"><p>Protein <span class="mono">{_esc(protein_id)}</span> not in DB.</p></div>',
                )
                return
            feats = p.get("features") or {}
            feat_html = "".join(
                f'<div><span class="k">{_esc(k)}</span><span class="v">{_esc(feats[k])}</span></div>'
                for k in sorted(feats.keys())
            )
            uniprot = f"https://www.uniprot.org/uniprotkb/{_esc(p['protein_id'])}/entry"
            go = ", ".join(_esc(g) for g in (p.get("go_terms_list") or [])[:24]) or "—"
            kw = ", ".join(_esc(g) for g in (p.get("keywords_list") or [])[:24]) or "—"
            body = f"""
            <div class="detail-hero">
              <div class="card">
                <h2 class="mono">{_esc(p['protein_id'])}</h2>
                <p><strong>{_esc(p.get('gene_name') or '—')}</strong>
                   · <span class="tier {_esc(p['tier'])}">{_esc(p['tier'])}</span>
                   · length {_esc(p.get('length'))}
                   · <a href="{uniprot}" rel="noopener" target="_blank">UniProt</a></p>
                <p>{_esc(p.get('protein_name') or 'No UniProt protein name in DB.')}</p>
                <p><span class="muted">GO:</span> {go}</p>
                <p><span class="muted">Keywords:</span> {kw}</p>
                <p class="muted mono">hash={_esc(p.get('sequence_hash'))} · source={_esc(p.get('source'))}</p>
              </div>
              <div class="card">
                <h3>Scores</h3>
                <p>PlantLLPS: <strong>{_fmt_score(p.get('score'))}</strong></p>
                <p>Companion: {_fmt_score(p.get('score_companion'))}</p>
                <p>Logistic: {_fmt_score(p.get('score_logistic'))}</p>
                <p>Percentile: {_fmt_pct(p.get('percentile'))}</p>
                <p>MAP65-tail-like: {_fmt_score(p.get('map65_tail_like'))}</p>
                <p>Disorder proxy: {_fmt_score(p.get('frac_disorder_proxy'))}</p>
                <p>metapredict frac: {_fmt_score(p.get('frac_disorder_metapredict'))}</p>
                <p>AF mean pLDDT: {_fmt_score(p.get('af_mean_plddt'))}
                   · frac&lt;50: {_fmt_score(p.get('af_frac_very_low'))}</p>
                <p>PLAAC-approx NLLR: {_fmt_score(p.get('plaac_nllr_approx'))}</p>
                <p>PSAP: {_fmt_score(p.get('score_psap'))}
                   · LLPhyScore %ile: {_fmt_score(p.get('score_llphyscore_percentile'))}</p>
                <p>Ath: <span class="mono">{_esc(p.get('ortholog_ath') or '—')}</span>
                   ({_esc(p.get('ortholog_ath_gene') or '')})
                   · Rice: <span class="mono">{_esc(p.get('ortholog_osa') or '—')}</span>
                   ({_esc(p.get('ortholog_osa_gene') or '')})
                   · family {_esc(p.get('cross_crop_family') or '—')}</p>
                <p>Expr. mean TPM: {
                    f"{float(p['expr_mean_tpm']):.2f}"
                    if p.get("expr_mean_tpm") is not None
                    else "—"
                }{
                    f' (<span class="mono">{_esc(p.get("expr_experiment"))}</span>)'
                    if p.get("expr_experiment")
                    else ""
                }</p>
                <p>Expr. max TPM: {
                    f"{float(p['expr_max_tpm']):.2f}"
                    if p.get("expr_max_tpm") is not None
                    else "—"
                }</p>
              </div>
            </div>
            {
                (
                    '<div class="card" style="margin-top:1rem"><h3>Expression Atlas TPM by condition</h3>'
                    '<div class="features">'
                    + "".join(
                        f'<div><span class="k">{_esc(k)}</span><span class="v">{_esc(v)}</span></div>'
                        for k, v in sorted((p.get("expr_tpm") or {}).items())
                    )
                    + "</div></div>"
                )
                if p.get("expr_tpm")
                else ""
            }
            <div class="card" style="margin-top:1rem">
              <h3>Features</h3>
              <div class="features">{feat_html}</div>
            </div>
            """
            self._html(200, p["protein_id"], body)

        def _score_form(self, result_html: str = "") -> None:
            body = f"""
            <div class="card">
              <h2>Score a sequence</h2>
              <p class="muted">Paste FASTA or raw amino-acid sequence. Uses the local PlantLLPS model
              (not written to the proteome DB unless you use the CLI).</p>
              <form method="post" action="/score">
                <p><label>Optional ID<br/>
                  <input type="text" name="protein_id" placeholder="QUERY1" style="min-width:12rem"/>
                </label></p>
                <p><label>Sequence<br/>
                  <textarea name="sequence" placeholder=">QUERY1&#10;MKT..."></textarea>
                </label></p>
                <button type="submit">Score</button>
              </form>
            </div>
            {result_html}
            """
            self._html(200, "Score FASTA", body)

        def _score_submit(self, form: Dict[str, List[str]]) -> None:
            raw = (form.get("sequence") or [""])[0]
            pid = (form.get("protein_id") or ["QUERY"])[0].strip() or "QUERY"
            seq = _extract_sequence(raw)
            if len(seq) < 20:
                return self._score_form(
                    '<div class="card err"><p>Need a protein sequence of at least 20 amino acids.</p></div>'
                )
            model = state.model()
            scored = model.score(seq)
            feat_html = "".join(
                f'<div><span class="k">{_esc(k)}</span><span class="v">{_esc(scored.get(k))}</span></div>'
                for k in sorted(k for k in scored.keys() if k not in {
                    "score_plantllps", "score_companion", "score_logistic", "tier"
                })
            )
            result = f"""
            <div class="card" style="margin-top:1rem">
              <h3>Result for <span class="mono">{_esc(pid)}</span></h3>
              <p>Length {len(seq)} · Tier
                 <span class="tier {_esc(scored.get('tier'))}">{_esc(scored.get('tier'))}</span></p>
              <p>PlantLLPS: <strong>{_fmt_score(scored.get('score_plantllps'))}</strong>
                 · companion {_fmt_score(scored.get('score_companion'))}
                 · logistic {_fmt_score(scored.get('score_logistic'))}</p>
              <div class="features">{feat_html}</div>
            </div>
            """
            self._score_form(result)

    return Handler


def _extract_sequence(raw: str) -> str:
    lines = []
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith(">"):
            continue
        lines.append(line)
    seq = "".join(lines) if lines else raw.strip()
    return re.sub(r"[^A-Za-z]", "", seq).upper()


def run_server(
    *,
    host: str = DEFAULT_HOST,
    port: int = DEFAULT_PORT,
    db_path: Path | str = DEFAULT_DB,
    model_path: Optional[Path | str] = None,
) -> None:
    db_path = Path(db_path)
    if not db_path.exists():
        raise FileNotFoundError(
            f"SQLite DB not found: {db_path}\n"
            "Run `python3 -m plantllps score-proteome` first, or pass --db."
        )
    state = AppState(db_path, Path(model_path) if model_path else None)
    handler = make_handler(state)
    httpd = ThreadingHTTPServer((host, port), handler)
    url = f"http://{host}:{port}/"
    preview = f"http://127.0.0.1:{port}/" if host in {"0.0.0.0", "::"} else url
    print(f"GrapeCondensatome UI → {preview}", flush=True)
    print(f"Bound {host}:{port}", flush=True)
    print(f"SQLite: {db_path}", flush=True)
    print("Ctrl+C to stop.", flush=True)
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print("\nStopped.", flush=True)
    finally:
        httpd.server_close()
