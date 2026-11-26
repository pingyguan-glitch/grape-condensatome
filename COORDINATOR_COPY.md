# Copy to coordinator `/agent/grape-condensatome`

From Project store:

```bash
rm -rf /agent/grape-condensatome
cp -a /cursor/stores/bc-fbfd25f5-3e47-43bf-8b41-3a65ef77fb51/internal/grape-condensatome-app /agent/grape-condensatome
cd /agent/grape-condensatome
npm install
npm start
# preview http://127.0.0.1:8765/
```

Then publish-app: `project_dir=grape-condensatome`, `preview_url=http://127.0.0.1:8765/`.
