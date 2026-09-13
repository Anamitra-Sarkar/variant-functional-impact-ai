#!/bin/sh
set -e
if [ "$MODEL_RELEASE_APPROVED" = "true" ] && [ -n "$APPROVED_ARTIFACT_REVISION" ] && [ -n "$MODEL_ARTIFACT_REPO_ID" ]; then
  echo "[entrypoint] Downloading approved model artifact from $MODEL_ARTIFACT_REPO_ID"
  python3 -c "
import os, shutil
from huggingface_hub import hf_hub_download
dest = os.environ.get('MODEL_ARTIFACT_PATH', 'artifacts/model.pkl')
os.makedirs(os.path.dirname(dest) or '.', exist_ok=True)
downloaded = hf_hub_download(
    repo_id=os.environ['MODEL_ARTIFACT_REPO_ID'],
    filename='model.pkl',
    repo_type='model',
    token=os.environ.get('HF_TOKEN'),
)
shutil.copy(downloaded, dest)
print('[entrypoint] Download complete ->', dest)
"
else
  echo "[entrypoint] No approved release configured; starting in fail-closed abstention mode."
fi
exec uvicorn backend.app:app --host 0.0.0.0 --port 8000
