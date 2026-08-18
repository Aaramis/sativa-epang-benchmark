# sativa-epang-benchmark

Does replacing SATIVA's placement engine with EPA-ng change its answers?

![Runtime and agreement](results/figures/fig_sativa_speedup.png)

Three results, measured against unmodified SATIVA v0.9.3 on ITS alignments of 100 to 5402
sequences:

- **Speed.** 257 s to 12.7 s at 1600 sequences on 8 threads, 1536 s to 87 s at 5402 on 16.
  The gain is parallelism over queries; RAxML gains little from threads on a 242 column
  alignment.
- **Agreement.** 0.96 recall and 0.94 precision at 800 sequences with K=25 folds, against
  0.92 / 0.90 when unmodified SATIVA merely changes its own substitution model.
- **Detection.** Of 33 mislabels injected into three clades, unmodified SATIVA finds 28 and
  this version 33. Every miss happens when at least two sequences share the wrong label,
  which a strict leave-one-out cannot see and a fold can.

**Results and interpretation: [`RESULTS.md`](RESULTS.md).** The code: [`sativa_epang/`](sativa_epang/),
also published as a branch on a fork of `amkozlov/sativa`. The rest of this file says what
is here and how to run it.

SATIVA is Kozlov et al. 2016; EPA-ng is Barbera et al. 2019. This work is a benchmark of a
modification, not a new method.

## Layout

```
config.sh              paths, interpreters, SATIVA parameters. Everything sources this.
sativa_epang/          the modified SATIVA the report describes, with PROVENANCE.md
env/                   the two unmodified upstream checkouts (see Setup)
data_ascii/            the size gradient, n = 100 to 5402
data_injected/         the same subsets with k sequences deliberately relabelled
scripts/               eight scripts, below
results/
  runs/<dataset>__<condition>__rep<i>/   one measured run: run.json, .mis, .log
  runs_injected/                         the same, for the injected mislabel test
  summary.tsv                            one row per condition and size
  concordance.tsv                        per sequence agreement
  concordance_details.tsv                one row per disagreeing sequence
  concordance_vs_gamma*.tsv              the same, against the GTRGAMMA pinned reference
  figures/fig_sativa_speedup.{png,pdf}   the figure in RESULTS.md
  report.html                            the report as a standalone page
```

`env/` and the compiled RAxML binaries are gitignored: both are rebuilt by the setup steps
below, and neither belongs in a copy of this folder.

## Scripts

| | |
|---|---|
| `01_make_datasets.py` | Nested random subsets of one alignment, written as SATIVA's `.fasta` + `.tax`. `--ascii-taxonomy` transliterates taxon names, which the 2017 python 2 SATIVA needs. |
| `02_bench.py` | One measured run of one condition. Owns the condition table: which SATIVA, which engine, how many threads, how many folds, which environment overrides. |
| `03_plot.py` | `summary.tsv` and the three panel figure. |
| `04_run_all.sh` | Every measurement in the report, in five phases. About 3 h on 32 cores. |
| `05_concordance.py` | Which sequences each version flags, against a reference run. `--details` writes one row per disagreement. |
| `06_cutoff_sensitivity.py` | The same agreement recomputed at stricter confidence cutoffs, from the `.mis` files already on disk. |
| `07_injected_mislabels.py` | `build` writes the datasets with known wrong labels, `score` reads the runs back. |
| `08_build_report.py` | Renders `RESULTS.md` into `results/report.html`, figure inlined, so the report can be shared on its own. |
| `_timeit.py` | Runs one command in its own process group and reports wall time and peak RSS. |

## Setup

The two unmodified references are git checkouts of `amkozlov/sativa`:

```bash
(cd sativa_epang/raxml && make)                                                 # RAxML 8.2.3
git clone https://github.com/amkozlov/sativa.git env/sativa_upstream          # v0.9.3, python 3
git -C env/sativa_upstream worktree add --detach ../sativa_upstream_py2 8c31962 # last python 2 commit
cp sativa_epang/raxml/raxmlHPC8-* env/sativa_upstream/raxml/
cp sativa_epang/raxml/raxmlHPC8-* env/sativa_upstream_py2/raxml/
mkdir -p env/sativa_upstream/tmp env/sativa_upstream_py2/tmp
conda create -y -n sativa_py2 -c conda-forge python=2.7
```

Copying the RAxML binaries matters: all three versions must call the same ones, so the
comparison never changes binary. The python 3 environment (`PY3` in `config.sh`) needs
`ete3` and `epa-ng` 0.3.8, for instance:

```bash
conda create -y -n sativa_epang -c conda-forge -c bioconda python=3.10 ete3 epa-ng=0.3.8
```

`config.sh` holds no site specific path. Put yours in `config.local.sh`, which it sources
if present and which `.gitignore` skips:

```bash
SATIVA_BENCH_PY3_ENV="/path/to/conda/envs/sativa_epang"      # python 3, ete3, epa-ng
SATIVA_BENCH_PY2_ENV="/path/to/conda/envs/sativa_py2"        # bare python 2.7
SATIVA_BENCH_PLOT_PY="/path/to/conda/envs/plot/bin/python"   # matplotlib, markdown, pillow
SATIVA_BENCH_SOURCE_MSA="/path/to/alignment.fasta"
```

## Running

```bash
source config.sh
"$PY3" scripts/01_make_datasets.py --source "$SOURCE_MSA" --outdir data_ascii \
       --sizes 100,200,400,800,1600,full --seed 42 --ascii-taxonomy
./scripts/04_run_all.sh
"$PY3" scripts/07_injected_mislabels.py build
"$PY3" scripts/02_bench.py --datasets-json data_injected/datasets.json \
       --runs-dir results/runs_injected --datasets all \
       --conditions upstream_py3_raxml_T1,port_py3_epang_T8,port_py3_epang_T8_K25 \
       --repeat 1 --untimed
"$PY3" scripts/07_injected_mislabels.py score
```

Two rules the runner enforces, both learned the hard way.

**Timed runs work on local disk.** The project lives on a network filesystem where 200 MB
take 14.0 s to write against 0.14 s locally. Measured there, the same run varied by 60 %;
measured on local disk, by 2 %. `02_bench.py` puts the working directory and SATIVA's own
`-tmpdir` under `--scratch` and copies back only `.mis` and `.log`.

**Timed runs do not share the machine.** Runs that answer an agreement question pass
`--untimed`, which records `timing_trusted: false` in `run.json`; `03_plot.py` refuses to
draw those seconds. Several such batches can run at once on disjoint cores with `--cpus`.

## Adding a condition

Conditions live in one table at the top of `02_bench.py`. A condition names an
implementation, an engine, a thread count, and optionally a fold count, environment
overrides, or a `sativa.cfg` override:

```python
"port_py3_epang_T8_K25": {"impl": "port", "engine": "epang", "threads": 8, "folds": 25},
"upstream_py3_raxml_T8_gamma": {"impl": "upstream_py3", "engine": "raxml", "threads": 8,
                                "cfg": {"raxml_model": "GTRGAMMA",
                                        "epa_use_heuristic": "FALSE"}},
```

The environment variables the EPA-ng version understands are listed in
[`sativa_epang/PROVENANCE.md`](sativa_epang/PROVENANCE.md).
