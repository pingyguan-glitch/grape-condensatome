# Public hosting notes

GrapeCondensatome is a **Node (Express) + SQLite** app (`server.js`, default port **8765**).
It is **not** a static site: GitHub Pages cannot serve the SQLite-backed API/UI.

## Current status (honest)

| Mode | Status |
|------|--------|
| **Local server** | Supported — `npm install && npm start` → `http://127.0.0.1:8765/` |
| **Confirmed public URL** | **None yet.** Do not invent or paste placeholder hostnames as live. |
| Cursor **publish-app** | **Proposed** (slug `grape-condensatome`, proposal awaiting owner confirm + terms). Not live until confirmed — then paste the issued URL here. |
| Railway / Render / Fly | Configs in-repo (`railway.toml`, `render.yaml`, `fly.toml`, `Dockerfile`). No platform tokens in the Cloud Agent; deploy requires your account login — see steps below. |

When a real public HTTPS URL exists, replace the table row above and the README badges in the same PR. Health check path: `GET /api/health`.

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
pack stays under the limit (SQLite + FASTA remain). If a pack still exceeds 50 MiB,
use Render/Railway/Fly instead.

From a Cloud Agent / Cursor environment with this repo committed **and pushed**:

1. Start the preview: `npm install && npm start` → `http://127.0.0.1:8765/`
2. Call Cursor **publish-app** on project dir `.` with that preview URL (suggested slug `grape-condensatome`).
3. Confirm in the Cursor publish dialog (visibility: **public** for NAR/database citation; private / team also available).
4. Only after confirm: paste the issued URL into this file’s status table and the README badges.

If publish-app is blocked (pack size, unpushed git, security review, or platform refusal), use a
container host below.

## Deploy configs (in this repo)

| File | Platform | Role |
|------|----------|------|
| `Dockerfile` | All | `node:22-bookworm-slim`; `HOST=0.0.0.0` `PORT=8080` |
| `render.yaml` | Render | Blueprint free Web Service; build `npm ci`; start `npm start`; health `/api/health` |
| `railway.toml` | Railway | Dockerfile build; health `/api/health` |
| `fly.toml` | Fly.io | Dockerfile; `force_https`; health check `/api/health` |
| `.dockerignore` | Docker | Keeps image lean (excludes `dist/`, caches) |

### Exact user steps — Render (recommended free path)

1. Log in at [https://dashboard.render.com](https://dashboard.render.com) (GitHub OAuth).
2. **New → Blueprint** → connect `pingyguan-glitch/grape-condensatome` → apply `render.yaml`.
   - Or **New → Web Service** → connect the repo → Build `npm ci --omit=dev` → Start `npm start` → Node **22** → Health `/api/health`.
3. Wait for the first deploy (SQLite ~100 MiB ships with the repo).
4. Open the issued `https://….onrender.com` URL; verify `…/api/health` returns `"ok": true`.
5. Paste that URL into this file + README badges; commit.

Free tier may spin down after idle (cold start tens of seconds).

### Exact user steps — Railway

1. Log in at [https://railway.app](https://railway.app) (GitHub).
2. **New Project → Deploy from GitHub** → select this repo (uses `railway.toml` / `Dockerfile`).
   - Or CLI: `npx --yes @railway/cli login && npx --yes @railway/cli init && npx --yes @railway/cli up`
3. Generate a public HTTPS domain in the service Settings.
4. Verify `/api/health`; paste URL into this file + README.

### Exact user steps — Fly.io

1. Install CLI: [https://fly.io/docs/hands-on/install-flyctl/](https://fly.io/docs/hands-on/install-flyctl/)
2. `fly auth login`
3. `fly apps create grape-condensatome` (rename in `fly.toml` if taken)
4. `fly deploy`
5. Verify `https://grape-condensatome.fly.dev/api/health`; paste URL into this file + README.

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

Until a public URL exists, README uses **pending** shield placeholders (not fake hosts).
After deploy, point shields at the live host, e.g.:

```markdown
[![Public atlas](https://img.shields.io/badge/atlas-LIVE_URL-0B6E4F)](https://YOUR-HOST/)
[![Health](https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2FYOUR-HOST%2Fapi%2Fhealth&query=%24.ok&label=health&prefix=&suffix=)](https://YOUR-HOST/api/health)
```

Optional: register the public URL with [UptimeRobot](https://uptimerobot.com/) (or similar) on `GET /api/health` and add their badge.

## Related

- Zenodo data deposit (weights + SQLite; **no invented DOI**): [`ZENODO.md`](./ZENODO.md)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
