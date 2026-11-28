# GrapeCondensatome — Node Express + SQLite atlas (not for GitHub Pages).
# Recommended host: Fly.io (see docs/DEPLOY.md). Also works on Render / Railway / any Docker host.
FROM node:22-bookworm-slim

WORKDIR /app

# Dependencies first for layer cache
COPY package.json package-lock.json ./
RUN npm ci --omit=dev

# App source + bundled data (model weights, SQLite, public/, examples)
COPY . .

# Fail the image build if production scoring/atlas assets are missing
RUN test -f data/model_formal.json \
  && test -f data/grape_condensatome.sqlite \
  && test -f public/style.css \
  && test -f public/examples/score_demo.fasta \
  && test -f VERSION \
  && node -e "JSON.parse(require('fs').readFileSync('data/model_formal.json','utf8'))"

ENV HOST=0.0.0.0 \
    PORT=8080 \
    NODE_ENV=production

EXPOSE 8080

# Liveness: Express must be bound (see server.js early listen)
HEALTHCHECK --interval=30s --timeout=5s --start-period=90s --retries=3 \
  CMD node -e "fetch('http://127.0.0.1:'+(process.env.PORT||8080)+'/api/health').then(r=>process.exit(r.ok?0:1)).catch(()=>process.exit(1))"

USER node
CMD ["node", "server.js"]
