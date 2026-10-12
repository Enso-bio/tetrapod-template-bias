"""Shared paths, cached downloads with provenance, and small helpers (reptile and amphibian study)."""
from pathlib import Path
import datetime as dt
import hashlib
import json
import subprocess
import threading
import time

import requests

ROOT = Path(__file__).resolve().parent.parent
DATA = ROOT / 'data'
RAW = DATA / 'raw'
RESULTS = ROOT / 'results'
FIGURES = ROOT / 'figures'
TOOLS = ROOT / 'tools'
SEED = 20261009
for folder in (DATA, RAW, RESULTS, FIGURES):
    folder.mkdir(parents=True, exist_ok=True)

USER_AGENT = 'ComparativeLongevityStudy/1.0 (academic, non-commercial)'
_lock = threading.Lock()
_provenance_path = DATA / 'provenance.jsonl'


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def record(url, path):
    entry = {'url': url, 'path': str(Path(path).relative_to(ROOT)).replace('\\', '/'),
             'bytes': Path(path).stat().st_size, 'sha256': sha256(path),
             'retrieved_utc': dt.datetime.now(dt.timezone.utc).isoformat(timespec='seconds')}
    with _lock, open(_provenance_path, 'a', encoding='utf-8') as f:
        f.write(json.dumps(entry) + '\n')


def fetch(url, path, headers=None, timeout=120, tries=6, ok_404=False, params=None):
    """Download url to path once; later calls read the cached file. Returns bytes or None on 400/404."""
    path = Path(path)
    if path.exists():
        return path.read_bytes()
    path.parent.mkdir(parents=True, exist_ok=True)
    h = {'User-Agent': USER_AGENT}
    h.update(headers or {})
    for attempt in range(tries):
        try:
            r = requests.get(url, headers=h, timeout=timeout, params=params)
            if r.status_code in (400, 404) and ok_404:
                return None
            if r.status_code == 429:
                time.sleep(float(r.headers.get('Retry-After', 2 * (attempt + 1))))
                continue
            r.raise_for_status()
            tmp = path.with_suffix(path.suffix + '.part')
            tmp.write_bytes(r.content)
            tmp.replace(path)
            record(r.url, path)
            return path.read_bytes()
        except requests.RequestException:
            if attempt == tries - 1:
                raise
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f'download failed after {tries} attempts (rate limited): {url}')


def fetch_json(url, path, **kw):
    b = fetch(url, path, headers={'Accept': 'application/json'}, **kw)
    return None if b is None else json.loads(b)


def curl(url, path, tries=5):
    """Large downloads with curl (resumable, streaming)."""
    path = Path(path)
    if path.exists():
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + '.part')
    for _ in range(tries):
        p = subprocess.run(['curl.exe', '-s', '-S', '-f', '-L', '-C', '-', '-A', USER_AGENT, '-o', str(tmp), url],
                           capture_output=True, text=True)
        if p.returncode == 0:
            tmp.replace(path)
            record(url, path)
            return path
        time.sleep(3)
    raise RuntimeError(f'download failed {url}: {p.stderr}')


def binomial(name):
    """'Chrysemys picta bellii' -> 'chrysemys_picta'."""
    return '_'.join(name.replace('_', ' ').lower().split()[:2])
