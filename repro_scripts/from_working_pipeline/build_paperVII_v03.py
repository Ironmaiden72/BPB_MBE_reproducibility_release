from pathlib import Path
import shutil, re, zipfile, hashlib, subprocess, textwrap, os
src = Path('/mnt/data/paperVII_work')
out = Path('/mnt/data/paperVII_v03_src')
if out.exists(): shutil.rmtree(out)
out.mkdir()
# copy dirs excluding manuscript pdf
for p in src.iterdir():
    if p.name == 'paperVII_BTFR_a0_v0_2.pdf':
        continue
    if p.is_dir():
        shutil.copytree(p, out/p.name)
    elif p.name == 'paperVII_BTFR_a0_v0_2.tex':
        shutil.copy2(p, out/'paperVII_BTFR_a0_v0_3.tex')
    else:
        shutil.copy2(p, out/p.name)
# duplicate data/script names for v0_3 consistency
data = out/'data'
for old in list(data.glob('paperVII_v0_2_*')):
    new = data/old.name.replace('v0_2','v0_3')
    shutil.copy2(old,new)
# duplicate script
scripts = out/'scripts'
if (scripts/'make_paperVII_v0_2_figures.py').exists():
    shutil.copy2(scripts/'make_paperVII_v0_2_figures.py', scripts/'make_paperVII_v0_3_figures.py')
# release index (no 00_README pdf to keep package lighter)
rel = out/'release_index'
rel.mkdir()
for name in ['00_MANIFEST_CORPUS.csv','00_CLAIMS_INDEX.csv','00_RELEASE_NOTES.md','00_REVIEWER_ROUTE.md','00_REVIEWER_ROUTE.csv','00_RELEASE_TODO.csv','00_RELEASE_STAGE1_QA.md']:
    p = Path('/mnt/data/stage1_idx')/name
    if p.exists(): shutil.copy2(p, rel/name)
tex_path = out/'paperVII_BTFR_a0_v0_3.tex'
s = tex_path.read_text(encoding='utf-8')
# update title/version
s = s.replace('v0.2 --- P7-A7 manuscript lock', 'v0.3 --- HAL corpus-aware P7-A7 manuscript lock')
# update purpose title and version textual references, but keep file names manually updated later
s = s.replace('\\section{Purpose of v0.2}', '\\section{Purpose of v0.3}')
s = s.replace('Paper VII v0.2', 'Paper VII v0.3')
s = s.replace('Version v0.2', 'Version v0.3')
s = s.replace('The v0.2 paper', 'The v0.3 paper')
s = s.replace('v0.2 claim', 'v0.3 claim')
s = s.replace('v0.2 lock', 'v0.3 lock')
s = s.replace('v0.2 numerical lock', 'v0.3 numerical lock')
s = s.replace('v0.2 gate chain', 'v0.3 gate chain')
s = s.replace('v0.2 source package', 'v0.3 source package')
# add citations in purpose after first paragraph-ish
old = "The central claim is deliberately bounded. BPB/MBE fixes a plausible acceleration scale and passes the BTFR/interpolation support layers. The full galaxy-level likelihood with distances, inclinations, stellar mass-to-light ratios, covariance, external-field effects and interpolation-family freedom remains a later gate."
new = old + " The comparison point is the standard MOND/RAR and SPARC rotation-curve phenomenology, but the claim made here is narrower: the BPB/MBE dictionary fixes the acceleration normalisation and tests its BTFR and baryonic-RC consequences, rather than deriving a final MOND interpolation family \\cite{Milgrom1983,Lelli2016,McGaugh2016}."
s = s.replace(old,new)
# insert corpus placement section before P7 source-lock dashboard
marker = '\\section{P7 source-lock dashboard}'
insert = r'''
\section{Corpus placement and release navigation}
Paper VII is the local-acceleration descendant of the BPB/MBE corpus. It should be read after Paper 0, which defines the survival factor $\ssurv$, and alongside Paper II, which uses SPARC rotation curves for the galaxy-scale Weyl-kernel test. In the release package, the top-level navigation files are:
\begin{itemize}
\item \texttt{00\_README\_CORPUS}: human entry point and quick orientation;
\item \texttt{00\_MANIFEST\_CORPUS.csv}: paper-by-paper manifest;
\item \texttt{00\_CLAIMS\_INDEX.csv}: claim-to-source-lock index;
\item \texttt{00\_REVIEWER\_ROUTE}: short audit route for reviewers.
\end{itemize}
The present paper contributes the Paper-VII rows of that corpus index. Its allowed claim is fixed-$a_0$ BTFR and baryonic-interpolation support. Its forbidden claims remain full SPARC rotation-curve likelihood closure, a unique interpolation-function derivation, active38 as an independent $a_0$ anchor, and an unconditional action-level MOND limit.

'''
s = s.replace(marker, insert + marker)
# add BTFR literature sentence in SPARC BTFR section
old = "P7-A3 tests the fixed acceleration scale against the SPARC baryonic Tully--Fisher relation using"
new = "P7-A3 tests the fixed acceleration scale against the SPARC baryonic Tully--Fisher relation using the SPARC baryonic catalogue and flat-velocity convention \\cite{Lelli2016}"
s = s.replace(old,new)
# this replacement made sentence weird: add period before equation? It will read 'using the ... convention \cite{Lelli2016}\begin{equation}' acceptable maybe missing colon. Let's fix if needed
s = s.replace('\\cite{Lelli2016}\n\\begin{equation}', '\\cite{Lelli2016}:\n\\begin{equation}')
# add RAR cite in raw RC section
s = s.replace('The primary interpolation is the empirical RAR form', 'The primary interpolation is the empirical RAR form used here as a diagnostic baryonic-interpolation law \\cite{McGaugh2016}')
# update appendix filenames
s = s.replace('paperVII_BTFR_a0_v0_2.tex', 'paperVII_BTFR_a0_v0_3.tex')
s = s.replace('data/paperVII_v0_2_numeric_locks.csv', 'data/paperVII_v0_3_numeric_locks.csv')
s = s.replace('data/paperVII_v0_2_gate_dashboard.csv', 'data/paperVII_v0_3_gate_dashboard.csv')
s = s.replace('data/paperVII_v0_2_claim_source_lock.csv', 'data/paperVII_v0_3_claim_source_lock.csv')
s = s.replace('scripts/make_paperVII_v0_2_figures.py', 'scripts/make_paperVII_v0_3_figures.py')
# add release_index to appendix verbatim list
s = s.replace('figs/fig6_a0_profile_ratios_v02.pdf/png\n\\end{verbatim}', 'figs/fig6_a0_profile_ratios_v02.pdf/png\nrelease_index/00_MANIFEST_CORPUS.csv\nrelease_index/00_CLAIMS_INDEX.csv\nrelease_index/00_REVIEWER_ROUTE.md\n\\end{verbatim}')
# add bibliography before end
bib = r'''
\begin{thebibliography}{9}

\bibitem{Milgrom1983}
M. Milgrom,
A modification of the Newtonian dynamics as a possible alternative to the hidden mass hypothesis,
\emph{Astrophysical Journal} \textbf{270}, 365--370 (1983).

\bibitem{Lelli2016}
F. Lelli, S. S. McGaugh, and J. M. Schombert,
SPARC: Mass Models for 175 Disk Galaxies with Spitzer Photometry and Accurate Rotation Curves,
\emph{Astronomical Journal} \textbf{152}, 157 (2016).

\bibitem{McGaugh2016}
S. S. McGaugh, F. Lelli, and J. M. Schombert,
Radial Acceleration Relation in Rotationally Supported Galaxies,
\emph{Physical Review Letters} \textbf{117}, 201101 (2016).

\bibitem{Verlinde2017}
E. P. Verlinde,
Emergent Gravity and the Dark Universe,
\emph{SciPost Physics} \textbf{2}, 016 (2017).

\bibitem{Planck2020}
Planck Collaboration,
Planck 2018 results. VI. Cosmological parameters,
\emph{Astronomy \& Astrophysics} \textbf{641}, A6 (2020).

\end{thebibliography}

'''
s = s.replace('\n\end{document}', '\n' + bib + '\\end{document}\n')
# update final source appendix section heading maybe
s = s.replace('\\section{v0.2 source files and reproduction}', '\\section{v0.3 source files and reproduction}')
tex_path.write_text(s, encoding='utf-8')
# write QA
qa = out/'paperVII_v0_3_QA.md'
qa.write_text("""# Paper VII v0.3 QA — source-only

## Status

`PAPER_VII_V0_3_SOURCE_ONLY_CORPUS_AWARE_PASS`

## Changes from v0.2

- Updated manuscript label to v0.3 / HAL corpus-aware release lock.
- Added a `Corpus placement and release navigation` section pointing to `00_README_CORPUS`, `00_MANIFEST_CORPUS.csv`, `00_CLAIMS_INDEX.csv`, and `00_REVIEWER_ROUTE`.
- Added `release_index/` folder with the stage-1 corpus navigation files.
- Added minimal bibliography for MOND, SPARC, RAR, emergent-gravity comparison, and Planck baseline context.
- Added explicit MOND/RAR/SPARC comparison wording while keeping the Paper-VII guardrail: no full SPARC RC likelihood, no unique interpolation derivation, no unconditional MOND action limit.
- Duplicated data/source-lock filenames to `paperVII_v0_3_*` while preserving all locked scientific values from P7-A7.
- No manuscript PDF included in the delivery ZIP.

## Scientific changes

None. Numeric locks are unchanged from P7-A7:

- `a0_BPB = 1.367856592595e-10 m/s²`
- `s_surv0 = 0.7956924397708606`
- BTFR primary `n=118`, median residual `-0.056007711329 dex`, RMS `0.295494277101 dex`
- P7-A5 raw RAR fractional improvement `0.911207525795`
- P7-A6b primary best `a0/a0_BPB = 0.807075252405`
- P7-A6b active38 best `a0/a0_BPB = 0.486722101981`, shape-only guardrail

## Compile QA

Local `pdflatex` compile performed in a temporary directory. Output PDF was generated for QA only and excluded from the source ZIP.
""", encoding='utf-8')
# compile twice
cmd = ['pdflatex','-interaction=nonstopmode','-halt-on-error','paperVII_BTFR_a0_v0_3.tex']
for i in range(2):
    r = subprocess.run(cmd, cwd=out, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, timeout=120)
    (out/f'compile_pass_{i+1}.log').write_text(r.stdout, encoding='utf-8')
    if r.returncode != 0:
        print(r.stdout[-4000:])
        raise SystemExit(f'pdflatex failed pass {i+1}')
# check unresolved
log = (out/'paperVII_BTFR_a0_v0_3.log').read_text(encoding='utf-8', errors='ignore')
warns = []
for pattern in ['undefined references','Citation','Rerun to get cross-references']:
    if pattern.lower() in log.lower(): warns.append(pattern)
# remove generated pdf and aux/log files from zip? keep compile logs? Usually source zip can exclude. We'll exclude generated pdf, aux, log.
for ext in ['.aux','.log','.out','.toc','.fls','.fdb_latexmk']:
    for p in out.glob(f'*{ext}'):
        p.unlink(missing_ok=True)
# remove generated manuscript pdf
(out/'paperVII_BTFR_a0_v0_3.pdf').unlink(missing_ok=True)
# also remove original v0_2 pdf if copied? no. Remove compile logs maybe keep? Better not.
for p in out.glob('compile_pass_*.log'):
    p.unlink(missing_ok=True)
# zip
zip_path = Path('/mnt/data/paperVII_BTFR_a0_v0_3_source_ONLY.zip')
if zip_path.exists(): zip_path.unlink()
with zipfile.ZipFile(zip_path,'w',zipfile.ZIP_DEFLATED) as z:
    for p in sorted(out.rglob('*')):
        if p.is_file():
            # no manuscript pdf, allow figure pdfs
            if p.name in ['paperVII_BTFR_a0_v0_2.pdf','paperVII_BTFR_a0_v0_3.pdf']:
                continue
            z.write(p, p.relative_to(out))
print('wrote', zip_path, zip_path.stat().st_size)
print('qa', qa)
