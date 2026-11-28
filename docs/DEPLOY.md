# Deploy GrapeCondensatome (independent HTTPS host)

**Recommended path: [Fly.io](https://fly.io)** — anonymous HTTPS on `:443`, Docker image with bundled model + SQLite, no Cursor apps-gate.

Alternates (`render.yaml`, `railway.toml`) remain in-repo for cheap fallbacks; Render’s free tier spins down after idle (cold starts), so it is a weaker NAR demo host than Fly with `min_machines_running = 1`.

Agents **cannot** create your cloud account or mint deploy tokens. You (or a teammate with billing access) must log in once and run the commands below.

## What ships in the image

| Asset | Path | Notes |
|-------|------|-------|
| Model weights | `data/model_formal.json` | Bundled in git / Docker (tiny) |
| Proteome atlas | `data/grape_condensatome.sqlite` | ~100 MiB; bundled |
| Static UI | `public/` | CSS, sample FASTA |
| Health | `GET /api/health` | JSON; `ready` becomes true after sidecars warm |
| Score demo | `/score?example=1` | No login |

No boot-time weight download is required for the default demo. Optional large caches (AlphaFold raw, OrthoDB dumps) are dockerignored and not needed for `/score` or the main atlas.

## Credentials you must provide

| Platform | What you create | Where |
|----------|-----------------|-------|
| **Fly.io (recommended)** | Fly account + payment method if required by Fly’s current free/hobby rules; then `fly auth login` (browser) or a **Fly API token** (`fly tokens create` / dashboard) for CI | [fly.io/dashboard](https://fly.io/dashboard) |
| Render (alternate) | Render account (GitHub OAuth) | [dashboard.render.com](https://dashboard.render.com) |
| Railway (alternate) | Railway account (GitHub) | [railway.app](https://railway.app) |

Do **not** paste tokens into chat. Store them in your password manager or Fly/Render secret stores.

## Fly.io — exact steps

1. Install the CLI: [https://fly.io/docs/hands-on/install-flyctl/](https://fly.io/docs/hands-on/install-flyctl/)
2. Authenticate on **your** machine:
   ```bash
   fly auth login
   ```
3. From the repo root (this branch / `main` after merge):
   ```bash
   # If the default name is taken, edit app = "…" in fly.toml first
   fly apps create grape-condensatome
   fly deploy
   ```
4. Verify (replace hostname if you renamed the app):
   ```bash
   curl -sS https://grape-condensatome.fly.dev/api/health | head
   curl -sS -o /dev/null -w "%{http_code}\n" https://grape-condensatome.fly.dev/score?example=1
   curl -sS -o /dev/null -w "%{http_code}\n" https://grape-condensatome.fly.dev/license
   curl -sS -o /dev/null -w "%{http_code}\n" https://grape-condensatome.fly.dev/privacy
   ```
5. Record the public URL in `docs/HOSTING.md` and the README badges in the same follow-up PR. **Do not invent a URL** before `fly deploy` succeeds.
6. Optional custom domain (NAR-friendly): Fly dashboard → Certificates → add your DNS CNAME.

### Fly config knobs (already in `fly.toml`)

- `force_https = true` — HTTPS :443 termination
- `min_machines_running = 1`, `auto_stop_machines = "off"` — warm machine for proposal day; set `min_machines_running = 0` and `auto_stop_machines = "stop"` later to save credits
- Health check: `GET /api/health`, 90 s grace
- VM: 1 GB RAM shared CPU (SQLite + fruit-crop sidecars)

### ≥5-year maintenance story (honest)

- Keep the Fly app billed/active; renew DNS if you attach a custom domain.
- Point NAR reviewers at the Fly (or custom) URL, **not** `*.apps.grok.bot` (Cursor gate).
- Source of truth remains this GitHub repo + MIT `LICENSE`; redeploy with `fly deploy` after merges.
- Zenodo DOI for data/software citation is separate (`docs/ZENODO.md`) and not required to keep HTTPS up.

## Render — alternate (Blueprint)

1. Log in at [https://dashboard.render.com](https://dashboard.render.com) (GitHub OAuth).
2. **New → Blueprint** → connect `pingyguan-glitch/grape-condensatome` → apply `render.yaml`.
3. Wait for the first deploy; open the issued `https://….onrender.com` URL.
4. Free plan cold-starts after idle — expect 30–60 s before `/api/health` answers.

## Railway — alternate

1. Log in at [https://railway.app](https://railway.app).
2. **New Project → Deploy from GitHub** (uses `railway.toml` / `Dockerfile`), or CLI `railway login && railway up`.
3. Generate a public HTTPS domain in service Settings.

## Local / Docker smoke test (no cloud account)

```bash
npm ci
HOST=0.0.0.0 PORT=8080 npm start
# → http://127.0.0.1:8080/api/health

docker build -t grape-condensatome .
docker run --rm -p 8080:8080 grape-condensatome
```

## Public pages to check after deploy

| Path | Expect |
|------|--------|
| `/` | Landing: MIT + privacy blurb; links to sample score / methods / limitations |
| `/score?example=1` | Sample FASTA scoring UI (no login) |
| `/methods` | Scoring methods |
| `/limitations` | Caveats |
| `/license` | Full MIT text |
| `/privacy` | FASTA retention / no public job listing |
| `/api/health` | `{"ok":true,...}` |

## Related

- Hosting status table: [`HOSTING.md`](./HOSTING.md)
- Web scoring caps: [`WEB_SCORING.md`](./WEB_SCORING.md)
