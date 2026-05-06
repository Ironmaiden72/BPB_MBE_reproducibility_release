#!/usr/bin/env python3
from pathlib import Path

p = Path('outputs/paperVI_EG_theory/paperVI_EG_T10_manifest.md')
if not p.exists():
    raise SystemExit(f'Missing {p}')

s = p.read_text(encoding='utf-8')
start = s.find('## Guardrail')
if start == -1:
    raise SystemExit('Could not find ## Guardrail section')

# Keep everything before Guardrail and replace the corrupted block to EOF.
prefix = s[:start].rstrip()
fixed = r'''

## Guardrail

This section states a projected-diagnostic closure:

```tex
r_{W\delta}=\sqrt{a_{\rm bg}}
```

It does not claim a completed action-level HR/BPB perturbation proof. The open field-equation target remains:

```tex
\frac{\langle P_m W,P_m W\rangle_{E_G}}{\langle W,W\rangle_{E_G}} = a_{\rm bg}
```
'''
p.write_text(prefix + fixed, encoding='utf-8')
print('PATCHED', p)
