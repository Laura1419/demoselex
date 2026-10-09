"""Extract the immutable peer-review PDF into physical-page-indexed derived text."""
from pathlib import Path
import hashlib
import json
from pypdf import PdfReader

source = Path('data/raw/akbari/41586_2026_10358_MOESM4_ESM.pdf')
expected = '8c94cddd5cd6c851403f68a4c3168dedfe5f5d29d1a973a808c0a0960365d108'
assert hashlib.sha256(source.read_bytes()).hexdigest() == expected
reader = PdfReader(source)
assert len(reader.pages) == 61
pages = [{'page': i + 1, 'text': page.extract_text()} for i, page in enumerate(reader.pages)]
assert all(page['text'].strip() for page in pages)
out = Path('data/derived/akbari_peer_review')
out.mkdir(parents=True, exist_ok=True)
(out / 'pages.json').write_text(json.dumps(pages, ensure_ascii=False, indent=2))
(out / 'text.txt').write_text('\n\n'.join(
    f"===== PDF PAGE {page['page']} =====\n{page['text']}" for page in pages))
print('PASS: verified source SHA256; extracted 61 nonempty pages with pypdf')
