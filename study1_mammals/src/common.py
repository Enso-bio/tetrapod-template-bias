"""Shared paths, cached downloads with provenance, and small helpers."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import threading
import time

import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
RAW = DATA / 'raw'
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
TOOLS = ROOT / 'tools'
SEED = 20261008
for folder in (DATA, RAW, RESULTS, FIGURES):
    folder.mkdir(parents=True, exist_ok=True)

USER_AGENT = 'MammalProteinStabilityStudy/1.0 (academic, non-commercial)'
_lock = threading.Lock()
_provenance_path = DATA / 'provenance.jsonl'


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _record(url, path):
    entry = {'url': url, 'path': str(Path(path).relative_to(ROOT)).replace('\\', '/'),
             'bytes': Path(path).stat().st_size, 'sha256': sha256(path),
             'retrieved_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}
    with _lock, open(_provenance_path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(entry) + '\n')


def fetch(url, path, headers=None, timeout=120, tries=6, ok_404=False):
    """Download url to path once; later calls read the cached file. Returns bytes or None on 404."""
    path = Path(path)
    if path.exists():
        return path.read_bytes()
    path.parent.mkdir(parents=True, exist_ok=True)
    h = {'User-Agent': USER_AGENT}
    h.update(headers or {})
    for attempt in range(tries):
        try:
            r = requests.get(url, headers=h, timeout=timeout)
            if r.status_code in (400, 404) and ok_404:  # Ensembl answers 400 for unknown symbols
                return None
            if r.status_code == 429:
                time.sleep(float(r.headers.get('Retry-After', 2 * (attempt + 1))))
                continue
            r.raise_for_status()
            tmp = path.with_suffix(path.suffix + '.part')
            tmp.write_bytes(r.content)
            tmp.replace(path)
            _record(url, path)
            return path.read_bytes()
        except requests.RequestException:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f'download failed after {tries} attempts (rate limited): {url}')


def fetch_json(url, path, **kw):
    b = fetch(url, path, headers={'Accept': 'application/json'}, **kw)
    return None if b is None else json.loads(b)


def binomial(name):
    """'canis_lupus_familiaris' -> 'canis_lupus'."""
    return '_'.join(name.lower().split('_')[:2])
