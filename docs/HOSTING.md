# Public hosting notes

GrapeCondensatome is a **Node (Express) + SQLite** app (`server.js`, default port **8765**).
It is **not** a static site: GitHub Pages cannot serve the SQLite-backed API/UI.

## Current status (honest)

| Mode | Status |
|------|--------|
| **Local server** | Supported — `npm install && npm start` → `http://127.0.0.1:8765/` |
| **Confirmed public URL** | **None.** Do not invent or paste placeholder URLs. |
| Cursor **publish-app** | May have been *proposed* in a Cloud Agent run; only live after the owner confirms the publish dialog + terms. Until then, treat hosting as **local-only**. |
| Railway / Render / Fly | Documented below; not deployed from this repo unless you deploy them. |

When a real public URL exists, record it here (and in the project store honesty table) in the same PR that confirms the deploy.

## Local preview (default)

```bash
npm install && npm start
# → http://127.0.0.1:8765/
# Browse /download or call:
#   curl 'http://127.0.0.1:8765/api/download?tier=A&format=tsv' -o tierA.tsv
#   curl 'http://127.0.0.1:8765/api/download?family=MAP65&format=fasta' -o map65.fasta
```

Runtime needs `data/UP000009183.fasta.gz` for FASTA downloads (indexed at startup).
Ship `data/grape_condensatome.sqlite` (~100 MiB at v1.1.0) with the app.

## Cursor publish-app

From a Cloud Agent / Cursor environment with this repo committed **and pushed**:

1. Start the preview: `npm install && npm start` → `http://127.0.0.1:8765/`
2. Call Cursor **publish-app** on project dir `.` with that preview URL.
3. Confirm in the Cursor publish dialog (visibility: private / team / public).
4. Only after confirm: paste the issued URL into this file’s status table.

If publish-app is blocked (unpushed git, security review, or platform refusal), use a
container host below.

## Alternatives (SQLite server)

| Option | Fit | Notes |
|--------|-----|-------|
| **Railway / Render / Fly.io** free tier | Good | Deploy Node 22+, ship `data/grape_condensatome.sqlite` (~100 MiB), set `PORT` |
| **Docker on a VPS** | Good | `node:22` image; bind `0.0.0.0:$PORT` |
| **GitHub Pages** | **No** | Static only; cannot run Express/SQLite |
| **Python `plantllps serve`** | Local / VPS | Stdlib HTTP UI; same SQLite |

### Railway / Render one-liners (after you create an account)

Railway (CLI; creates a project from the current directory):

```bash
npm install && npx --yes @railway/cli login && npx --yes @railway/cli init && npx --yes @railway/cli up
# Set start command to: npm start
# Ensure PORT is provided by the platform (server.js already reads process.env.PORT)
```

Render (GitHub-connected Web Service — dashboard is the usual path):

```bash
# In Render: New → Web Service → connect this GitHub repo
# Build:  npm install
# Start:  npm start
# Env:    NODE_VERSION=22   (or choose Node 22 runtime)
```

### Fly.io sketch

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

Egress: downloading UniProt / Expression Atlas caches at build time needs outbound HTTPS;
runtime browse mode is local SQLite only.

## Related

- Zenodo data deposit (weights + SQLite; **no invented DOI**): [`ZENODO.md`](./ZENODO.md)
- External tool pins: [`EXTERNAL_TOOLS.md`](./EXTERNAL_TOOLS.md)
