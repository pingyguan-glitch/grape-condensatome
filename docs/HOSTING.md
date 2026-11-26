# Public hosting notes

GrapeCondensatome is a **Node (Express) + SQLite** app (`server.js`, default port **8765**).
It is **not** a static site: GitHub Pages cannot serve the SQLite-backed API/UI.

## Cursor publish-app (preferred for this repo run)

From a Cloud Agent / Cursor environment with this repo committed **and pushed**:

1. Start the preview: `npm install && npm start` → `http://127.0.0.1:8765/`
2. Call Cursor **publish-app** on project dir `.` with that preview URL.
3. Confirm in the Cursor publish dialog (visibility: private / team / public).

If publish-app is blocked (unpushed git, security review, or platform refusal), use a
container host below. This document records the honest status after each attempt.

## Alternatives (SQLite server)

| Option | Fit | Notes |
|--------|-----|-------|
| **Railway / Render / Fly.io** free tier | Good | Deploy Node 22+, ship `data/grape_condensatome.sqlite` (~36 MB), set `PORT` |
| **Docker on a VPS** | Good | `node:22` image; bind `0.0.0.0:$PORT` |
| **GitHub Pages** | **No** | Static only; cannot run Express/SQLite |
| **Python `plantllps serve`** | Local / VPS | Stdlib HTTP UI; same SQLite |

Example Fly.io sketch:

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

## Local preview

```bash
npm install && npm start
# → http://127.0.0.1:8765/
```
