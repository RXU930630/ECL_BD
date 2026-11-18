"""Automatic, content-addressed JSON caches; never deserialize executable objects."""
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy


def canonical(value):
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(',', ':'), allow_nan=False)


def method_fingerprint():
    root = Path(__file__).parent
    digest = hashlib.sha256()
    for path in sorted(root.glob('*.py')):
        digest.update(path.name.encode()); digest.update(path.read_bytes())
    return {'source_sha256': digest.hexdigest(), 'numpy': np.__version__, 'scipy': scipy.__version__,
            'cache_format': 1}


class TableCache:
    def __init__(self, base, conditions):
        self.root = Path(base) / '.cache' / 'tables'
        self.conditions = {**conditions, 'method': method_fingerprint()}
        self.description = {}

    def entry(self, kind, extra=None):
        conditions = {**self.conditions, **(extra or {})}
        key = hashlib.sha256(canonical(conditions).encode()).hexdigest()
        return self.root / kind / key, conditions

    def load(self, kind, extra=None):
        folder, conditions = self.entry(kind, extra)
        try:
            data = json.loads((folder / 'table.json').read_text(encoding='utf-8'))
            if data['conditions'] != conditions:
                return None
            if hashlib.sha256(canonical(data['table']).encode()).hexdigest() != data['table_sha256']:
                return None
            return data['table']
        except (OSError, ValueError, KeyError, TypeError):
            return None

    def save(self, kind, table, extra=None):
        folder, conditions = self.entry(kind, extra)
        envelope = {'description': self.description, 'conditions': conditions, 'created_utc': datetime.now(timezone.utc).isoformat(),
                    'table_sha256': hashlib.sha256(canonical(table).encode()).hexdigest(), 'table': table}
        try:
            folder.mkdir(parents=True, exist_ok=True)
            # Atomic single-file publication: interrupted writes are never cache hits.
            with tempfile.NamedTemporaryFile(mode='w', encoding='utf-8', dir=folder, delete=False) as f:
                temporary = f.name
                json.dump(envelope, f, ensure_ascii=False, indent=2, allow_nan=False)
            os.replace(temporary, folder / 'table.json')
            summary=dict(conditions)
            for key in ('pressure_grid_pa','gap_grid_m'):
                if key in summary:
                    values=summary[key]
                    summary[key]={'minimum':min(values),'maximum':max(values),'points':len(values)}
            (folder/'CONDITIONS.txt').write_text(
                '缓存类型 / Table type: '+kind+'\n生成时间 / UTC: '+envelope['created_utc']+
                '\nPSD 文件名及单位（仅作说明）:\n'+json.dumps(self.description,ensure_ascii=False,indent=2)+
                '\n计算条件（SI）:\n'+json.dumps(summary,ensure_ascii=False,indent=2)+'\n',encoding='utf-8')
        except OSError as exc:
            print(f'缓存未保存（计算结果仍有效）：{exc}', flush=True)

    def describe(self, kind, hit, extra=None):
        folder, _ = self.entry(kind, extra)
        print(f'{kind}：'+('读取缓存' if hit else '计算并缓存（无完整匹配缓存）')+f'；条件与数据：{folder / "table.json"}', flush=True)
