# Public hosting notes

GrapeCondensatome is a **Node (Express) + SQLite** app (`server.js`, default port **8765**).
It is **not** a static site: GitHub Pages cannot serve the SQLite-backed API/UI.

## Recommended independent host

| Mode | Status |
|------|--------|
| **Recommended for NAR** | **Fly.io** — see exact click/CLI steps in [`DEPLOY.md`](./DEPLOY.md) |
| **Confirmed Fly public URL** | **Not yet** — requires the repo owner to `fly auth login` + `fly deploy` (no token in this environment) |
| Cursor publish-app | `https://grape-condensatome.apps.grok.bot` — **Cursor apps-gate**; not anonymous-demo ready |
| Local server | `npm install && npm start` → `http://127.0.0.1:8765/` |
| Render / Railway | In-repo alternates (`render.yaml`, `railway.toml`); Render free spins down after idle |

Do **not** invent hostnames. After a successful Fly (or other) deploy, paste the real `https://…` URL into this table and the README badges in the same PR.

### Why not Cursor `apps.grok.bot` for NAR?

NAR Web Server expects anonymous HTTPS on `:443` with a multi-year maintenance story. The Cursor apps gate shows an “Open this app” interstitial and is unsuitable as the public proposal-day demo. Use Fly (or another non-gated container host).

### Security review notes (follow-up)

Cursor publish security review flagged **resource_exhaustion** and **data_exposure** class findings for this public atlas (broad download/search surfaces over a full proteome SQLite). **Not fixed in this hosting pass** — track for a later hardening PR (rate limits, download caps, tighter static/data exposure). The atlas data are intentionally public research tables; secrets are not embedded. Submitted FASTA privacy: `/privacy` in the UI (request-memory only; no public job listing).

## Local preview (default)

```bash
npm install && npm start
# → http://127.0.0.1:8765/
# Browse /download or call:
#   curl 'http://127.0.0.1:8765/api/download?tier=A&format=tsv' -o tierA.tsv
#   curl 'http://127.0.0.1:8765/api/download?family=MAP65&format=fasta' -o map65.fasta
#   curl -s http://127.0.0.1:8765/api/health | head
```

Runtime needs `data/UP000009183.fasta.gz` for FASTA downloads (indexed at startup).
Ship `data/grape_condensatome.sqlite` (~100 MiB at v1.1.x) and `data/model_formal.json` with the app (both in the Docker image).

## Deploy configs (in this repo)

| File | Platform | Role |
|------|----------|------|
| [`DEPLOY.md`](./DEPLOY.md) | — | **Start here** — credentials + Fly steps |
| `Dockerfile` | All | `node:22-bookworm-slim`; `HOST=0.0.0.0` `PORT=8080`; HEALTHCHECK; non-root |
| `fly.toml` | **Fly.io (recommended)** | `force_https`; warm machine; health `/api/health` |
| `render.yaml` | Render | Blueprint free Web Service; health `/api/health` |
| `railway.toml` | Railway | Dockerfile build; health `/api/health` |
| `Procfile` | Render/Heroku-style | `web: npm start` |
| `.dockerignore` | Docker | Keeps image lean (excludes `dist/`, caches) |

## Cursor publish-app (legacy / gated)

Pack limit is **50 MiB gzipped**. Flat proteome TSV dumps are gitignored so the
pack stays under the limit (SQLite + FASTA remain).

Live gated instance: **https://grape-condensatome.apps.grok.bot**

Prefer Fly for proposal day. Only republish to Cursor if you still want an internal preview.

## Uptime / README badges

Until Fly is deployed, README may still mention the Cursor host for historical continuity — treat it as **gated**, not NAR-ready. After Fly goes live, point public-atlas / uptime shields at the Fly (or custom) URL and register it with [UptimeRobot](https://uptimerobot.com/) (or similar) on `/api/health`.

## Related

- Deploy steps + credentials: [`DEPLOY.md`](./DEPLOY.md)
- Zenodo data deposit (weights + SQLite; **no invented DOI**): [`ZENODO.md`](./ZENODO.md)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
- Web scoring: [`WEB_SCORING.md`](./WEB_SCORING.md)
