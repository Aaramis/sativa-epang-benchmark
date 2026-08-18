# Publishing this folder

Everything below is prepared. Nothing is public until you push.

## The two repositories

**`sativa-epang-benchmark`** (this folder). Already a git repository with one commit.

```bash
git remote add origin git@github.com:<you>/sativa-epang-benchmark.git
git push -u origin main
```

**The fork**, prepared in `../sativa-epang`. Fork `amkozlov/sativa` on GitHub first, then:

```bash
cd ../sativa-epang
git remote add fork git@github.com:<you>/sativa.git
git push fork epa-ng
```

The branch holds one commit on top of v0.9.3 (`68284c2`), so GitHub shows the whole
modification as a single readable diff: two hunks in `sativa.py`, one in `epac/config.py`,
one new `epac/epang_l1o.py`, one config line. `CHANGES-epa-ng.md` at the root of the branch
explains it and links back to this benchmark.

If your GitHub handle is not `agardette`, fix the two links: `CHANGES-epa-ng.md` in the
fork, and the "the code" line in this folder's `README.md`.

## Check before making them public

**Sequence redistribution.** `data_ascii/` and `data_injected/` hold 5 MB of ITS sequences
drawn from your curated database. If redistributing them is not clear cut, drop them and
keep the pipeline reproducible from a public source:

```bash
git rm -r --cached data_ascii data_injected
printf 'data_ascii/\ndata_injected/\n' >> .gitignore
```

`scripts/01_make_datasets.py` rebuilds both from any alignment, and `RESULTS.md` says which
one was used.

**Licence.** `LICENSE` is GPL 3, which SATIVA requires: `sativa_epang/` is a modified copy
of GPL 3 code, and `sativa_epang/PROVENANCE.md` states the modifications, as the licence
asks.

**Nothing site specific.** `config.local.sh` (your paths) and `env/` (two clones of a public
repository) are gitignored. The recorded runs were rewritten by
`scripts/10_sanitise_paths.py`, which replaces the interpreter, the benchmark root and the
scratch directory with `$PY3`, `./` and `$SCRATCH`. Re-run it after any new batch:

```bash
source config.sh
"$PY3" scripts/10_sanitise_paths.py --py3-env "$PY3_ENV" \
       --also data_ascii/datasets.json data_injected/datasets.json
git grep -l "/home/\|/pools/"   # should print nothing
```

## Repository description and topics

Description: *Does replacing SATIVA's placement engine with EPA-ng change its answers?
Speed, per-sequence agreement, and a ground-truth mislabel test.*

Topics: `sativa`, `epa-ng`, `phylogenetic-placement`, `taxonomy`, `benchmark`,
`reference-database`.
