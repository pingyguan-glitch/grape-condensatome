# Public hosting notes

GrapeCondensatome is a **Node (Express) + SQLite** app (`server.js`, default port **8765**).
It is **not** a static site: GitHub Pages cannot serve the SQLite-backed API/UI.

## Current status (honest)

| Mode | Status |
|------|--------|
| **Confirmed public URL** | **https://grape-condensatome.apps.grok.bot** (Cursor publish-app, public; confirmed 2026-11-27) |
| Health | Prefer `GET /api/health` after the app shell loads; first hit may show the Cursor apps-gate “Open this app” page |
| **Local server** | Still supported — `npm install && npm start` → `http://127.0.0.1:8765/` |
| Railway / Render / Fly | Configs in-repo as fallback / alternate hosts (`railway.toml`, `render.yaml`, `fly.toml`, `Dockerfile`) |

Do **not** invent other hostnames. If the public URL moves, update this table and the README badges in the same PR.

### Security review notes (follow-up)

Cursor publish security review flagged **resource_exhaustion** and **data_exposure** class findings for this public atlas (broad download/search surfaces over a full proteome SQLite). **Not fixed in this hosting pass** — track for a later hardening PR (rate limits, download caps, tighter static/data exposure). The atlas data are intentionally public research tables; secrets are not embedded.

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
Ship `data/grape_condensatome.sqlite` (~100 MiB at v1.1.x) with the app.

## Cursor publish-app

Pack limit is **50 MiB gzipped**. Flat proteome TSV dumps are gitignored so the
pack stays under the limit (SQLite + FASTA remain).

Live instance (this repo): **https://grape-condensatome.apps.grok.bot**

To republish after code changes: commit + push, keep a local preview on `:8765`, call **publish-app**, confirm if the platform asks.

If publish-app is blocked (pack size, unpushed git, security review hard-block, or platform refusal), use a
container host below.

## Deploy configs (in this repo)

| File | Platform | Role |
|------|----------|------|
| `Dockerfile` | All | `node:22-bookworm-slim`; `HOST=0.0.0.0` `PORT=8080` |
| `render.yaml` | Render | Blueprint free Web Service; build `npm ci`; start `npm start`; health `/api/health` |
| `railway.toml` | Railway | Dockerfile build; health `/api/health` |
| `fly.toml` | Fly.io | Dockerfile; `force_https`; health check `/api/health` |
| `.dockerignore` | Docker | Keeps image lean (excludes `dist/`, caches) |

### Exact user steps — Render (alternate)

1. Log in at [https://dashboard.render.com](https://dashboard.render.com) (GitHub OAuth).
2. **New → Blueprint** → connect `pingyguan-glitch/grape-condensatome` → apply `render.yaml`.
   - Or **New → Web Service** → connect the repo → Build `npm ci --omit=dev` → Start `npm start` → Node **22** → Health `/api/health`.
3. Wait for the first deploy (SQLite ~100 MiB ships with the repo).
4. Open the issued `https://….onrender.com` URL; verify health.
5. If this becomes the primary host, replace the URL in this file + README badges.

### Exact user steps — Railway (alternate)

1. Log in at [https://railway.app](https://railway.app) (GitHub).
2. **New Project → Deploy from GitHub** → select this repo (uses `railway.toml` / `Dockerfile`).
   - Or CLI: `npx --yes @railway/cli login && npx --yes @railway/cli init && npx --yes @railway/cli up`
3. Generate a public HTTPS domain in the service Settings.
4. Verify health; update this file + README if primary.

### Exact user steps — Fly.io (alternate)

1. Install CLI: [https://fly.io/docs/hands-on/install-flyctl/](https://fly.io/docs/hands-on/install-flyctl/)
2. `fly auth login`
3. `fly apps create grape-condensatome` (rename in `fly.toml` if taken)
4. `fly deploy`
5. Verify `https://grape-condensatome.fly.dev/api/health`; update docs if primary.

### Docker (any VPS)

```dockerfile
FROM node:22-bookworm-slim
WORKDIR /app
COPY package.json package-lock.json ./
RUN npm ci --omit=dev
COPY . .
ENV HOST=0.0.0.0 PORT=8080
EXPOSE 8080
CMD ["npm","start"]
```

```bash
docker build -t grape-condensatome .
docker run --rm -p 8080:8080 grape-condensatome
```

Egress: downloading UniProt / Expression Atlas caches at build time needs outbound HTTPS;
runtime browse mode is local SQLite only.

## Uptime / README badges

README points the public-atlas / uptime shields at **https://grape-condensatome.apps.grok.bot**.
Optional: register that URL with [UptimeRobot](https://uptimerobot.com/) (or similar). Note the Cursor apps-gate may return an HTML interstitial before the Express app answers `/api/health`, so prefer a simple HTTP-200 uptime check on `/` unless you automate the gate.

## Related

- Zenodo data deposit (weights + SQLite; **no invented DOI**): [`ZENODO.md`](./ZENODO.md)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
