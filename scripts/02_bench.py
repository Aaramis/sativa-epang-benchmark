#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Run the SATIVA speed benchmark: (dataset x condition) -> one run.json each.

Conditions differ in exactly one thing at a time:

  upstream_py2_raxml_T*  amkozlov/sativa @8c31962 (last python-2 commit), l1o = RAxML -f O
  upstream_py3_raxml_T*  amkozlov/sativa @68284c2 = v0.9.3, upstream's own python-3 port, RAxML -f O
  port_py3_raxml_T*      this repo's port, python3, l1o forced back to RAxML (SATIVA_L1O_ENGINE=raxml)
  port_py3_epang_T*      this repo's port, python3, l1o = EPA-ng k-fold, K=5 (the default)
  port_py3_epang_T8_K*   same, with K folds; K=N ("Kexact") prunes one single sequence per
                         run, which is the strict per-sequence leave-one-out SATIVA does --
                         the k-fold is the only approximation in the modified version, so
                         this sweep measures what exactness costs

upstream vs port at the same engine/threads isolates the python3 port;
raxml vs epang isolates the placement engine; T1 vs T8/T16 isolates threading.
Both implementations call the *same* RAxML 8.2.3 binaries for the reference
tree, so the shared `reftree` phase is comparable and can be subtracted.

Runs are strictly sequential -- a timing benchmark cannot share the machine
with itself.

Every run happens on LOCAL disk (``--scratch``, default /tmp/...): the project lives on
/pools, a network filesystem that writes 200 MB in ~14 s against ~0.14 s locally, so
SATIVA's own temp traffic (RAxML working files, .jplace) would dominate the measurement and
add several-fold run-to-run noise. Only the small artefacts (.mis, .log) are copied back
to results/runs/.
"""

import argparse
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
BENCH_ROOT = HERE.parent

CONDITIONS = {
    "upstream_py2_raxml_T1":  {"impl": "upstream_py2", "engine": "raxml", "threads": 1},
    "upstream_py2_raxml_T8":  {"impl": "upstream_py2", "engine": "raxml", "threads": 8},
    "upstream_py3_raxml_T1":  {"impl": "upstream_py3", "engine": "raxml", "threads": 1},
    "upstream_py3_raxml_T8":  {"impl": "upstream_py3", "engine": "raxml", "threads": 8},
    "port_py3_raxml_T1":      {"impl": "port",     "engine": "raxml", "threads": 1},
    "port_py3_raxml_T8":      {"impl": "port",     "engine": "raxml", "threads": 8},
    "port_py3_epang_T1":      {"impl": "port",     "engine": "epang", "threads": 1, "folds": 5},
    "port_py3_epang_T2":      {"impl": "port",     "engine": "epang", "threads": 2, "folds": 5},
    "port_py3_epang_T4":      {"impl": "port",     "engine": "epang", "threads": 4, "folds": 5},
    "port_py3_epang_T8":      {"impl": "port",     "engine": "epang", "threads": 8, "folds": 5},
    "port_py3_epang_T16":     {"impl": "port",     "engine": "epang", "threads": 16, "folds": 5},
    # Lining EPA-ng up with RAxML's placement rules, one at a time. SATIVA runs
    # `raxmlHPC -f O` with no preplacement heuristic (below 1000 taxa), RAxML-style branch
    # length optimisation and an accumulated-LWR cutoff of 0.999; EPA-ng defaults to a
    # two-phase heuristic, a sliding optimisation and (here) 0.99999. Each variable isolates
    # one of those, and *_raxmlmatch sets all of them plus the strict leave-one-out.
    "port_py3_epang_T8_noheur": {"impl": "port", "engine": "epang", "threads": 8,
                                 "env": {"SATIVA_EPANG_HEUR": "off"}},
    "port_py3_epang_T8_blo":    {"impl": "port", "engine": "epang", "threads": 8,
                                 "env": {"SATIVA_EPANG_BLO": "raxml"}},
    "port_py3_epang_T8_acc999": {"impl": "port", "engine": "epang", "threads": 8,
                                 "env": {"SATIVA_EPANG_ACC_LWR": "0.999"}},
    "port_py3_epang_T8_raxmlmatch": {"impl": "port", "engine": "epang", "threads": 8,
                                     "folds": 10 ** 6,
                                     "env": {"SATIVA_EPANG_HEUR": "off",
                                             "SATIVA_EPANG_BLO": "raxml",
                                             "SATIVA_EPANG_ACC_LWR": "0.999"}},
    # Above 500 taxa SATIVA drops its own RAxML to GTRCAT, and above 1000 it places
    # thoroughly on a fraction of the branches only -- so past those thresholds the
    # reference is a cheaper approximation and EPA-ng (always GTRGAMMA, every branch) is
    # not answering the same question. These conditions pin the reference to GTRGAMMA with
    # the heuristic off, which is what EPA-ng does, at the price of a much longer run.
    "upstream_py3_raxml_T8_gamma": {"impl": "upstream_py3", "engine": "raxml", "threads": 8,
                                    "cfg": {"raxml_model": "GTRGAMMA",
                                            "epa_use_heuristic": "FALSE"}},
    "upstream_py3_raxml_T1_gamma": {"impl": "upstream_py3", "engine": "raxml", "threads": 1,
                                    "cfg": {"raxml_model": "GTRGAMMA",
                                            "epa_use_heuristic": "FALSE"}},
    # sativa_epang/sativa.cfg now builds the reference tree under GTRGAMMA, so EPA-ng reads
    # a fitted alpha from RAxML_info. These two conditions restore the previous behaviour,
    # where a tree above 500 taxa was built with GTRCAT and its RAxML_info carried the
    # placeholder "alpha: 1.000000": they are what the before/after numbers in RESULTS.md
    # compare against.
    "port_py3_epang_T8_catinfo": {"impl": "port", "engine": "epang", "threads": 8,
                                  "cfg": {"raxml_model": "GTRCAT"}},
    "port_py3_epang_T8_K25_catinfo": {"impl": "port", "engine": "epang", "threads": 8,
                                      "folds": 25, "cfg": {"raxml_model": "GTRCAT"}},
    # Fold sweep: K=5 is the shipped default, K=N is the strict leave-one-out.
    "port_py3_epang_T8_K10":    {"impl": "port", "engine": "epang", "threads": 8, "folds": 10},
    "port_py3_epang_T1_K25":    {"impl": "port", "engine": "epang", "threads": 1, "folds": 25},
    "port_py3_epang_T8_K25":    {"impl": "port", "engine": "epang", "threads": 8, "folds": 25},
    "port_py3_epang_T8_K50":    {"impl": "port", "engine": "epang", "threads": 8, "folds": 50},
    "port_py3_epang_T8_K100":   {"impl": "port", "engine": "epang", "threads": 8, "folds": 100},
    "port_py3_epang_T8_Kexact": {"impl": "port", "engine": "epang", "threads": 8, "folds": 10 ** 6},
}

TIMING_RE = re.compile(
    r"elapsed time:\s*([\d.]+)\s*seconds\s*\(([\d.]+)s reftree,\s*([\d.]+)s leave-one-out\)"
)


def env_path(name, default=None):
    value = os.environ.get(name, default)
    if not value:
        sys.exit(f"missing environment variable {name} (source config.sh first)")
    return value


def count_mislabels(mis_path):
    """Number of flagged sequences in a SATIVA .mis file (';' lines are the header)."""
    if not mis_path.exists():
        return None
    n = 0
    with mis_path.open() as handle:
        for line in handle:
            if line.strip() and not line.startswith(";"):
                n += 1
    return n


def parse_sativa_log(log_path):
    """Pull SATIVA's own phase timings out of its log."""
    if not log_path.exists():
        return {}
    match = TIMING_RE.search(log_path.read_text(errors="replace"))
    if not match:
        return {}
    total, reftree, l1out = (float(g) for g in match.groups())
    return {"sativa_total_sec": total, "reftree_sec": reftree, "l1out_sec": l1out}


def build_command(cfg, condition, fasta, tax, run_name):
    spec = CONDITIONS[condition]
    interpreter, script = {
        "upstream_py2": (cfg["py2"], cfg["sativa_upstream_py2"]),
        "upstream_py3": (cfg["py3"], cfg["sativa_upstream_py3"]),
        "port":         (cfg["py3"], cfg["sativa_port"]),
    }[spec["impl"]]
    return [
        interpreter, str(script),
        "-s", fasta.name,
        "-t", tax.name,
        "-m", cfg["mode"],
        "-C", str(cfg["conf"]),
        "-x", cfg["taxcode"],
        "-T", str(spec["threads"]),
        "-p", str(cfg["seed"]),
        "-n", run_name,
        "-S",
    ]


def run_one(cfg, dataset, condition, timeout, force=False, rep=1):
    spec = CONDITIONS[condition]
    out_dir = Path(cfg["runs_dir"]) / f"{dataset['name']}__{condition}__rep{rep}"
    result_path = out_dir / "run.json"
    if result_path.exists() and not force:
        print(f"  [skip] {dataset['name']} / {condition} (already done)")
        return json.loads(result_path.read_text())

    run_name = f"{dataset['name']}_{condition}_rep{rep}"
    # Working directory and SATIVA's own temp directory both on local disk
    # (see the module docstring): on /pools the filesystem, not the algorithm, is
    # what gets timed.
    run_dir = Path(cfg["scratch"]) / "work" / run_name
    temp_dir = Path(cfg["scratch"]) / "tmp" / run_name
    for path in (run_dir, temp_dir):
        if path.exists():
            shutil.rmtree(path)
        path.mkdir(parents=True)

    fasta = run_dir / Path(dataset["fasta"]).name
    tax = run_dir / Path(dataset["tax"]).name
    shutil.copy(dataset["fasta"], fasta)
    shutil.copy(dataset["tax"], tax)

    command = build_command(cfg, condition, fasta, tax, run_name) + ["-tmpdir", str(temp_dir)]

    if spec.get("cfg"):
        # SATIVA reads its RAxML settings from sativa.cfg next to sativa.py; -c points it
        # at another one. raxml_home is resolved relative to the config file, so it has to
        # be spelled out absolutely here.
        script_dir = Path(build_command(cfg, condition, fasta, tax, run_name)[1]).resolve().parent
        cfg_path = run_dir / "sativa_override.cfg"
        lines = ["[raxml]",
                 f"raxml_home={script_dir / 'raxml'}",
                 "raxml_exec=run_raxml.sh",
                 "epa_load_optmod=true"]
        lines += [f"{key}={value}" for key, value in spec["cfg"].items()]
        cfg_path.write_text("\n".join(lines) + "\n")
        command += ["-c", str(cfg_path)]

    env = dict(os.environ)
    # epa-ng and the RAxML wrapper are resolved from PATH by SATIVA
    env["PATH"] = f"{cfg['py3_env']}/bin:" + env.get("PATH", "")
    if spec["impl"] == "port":
        env["SATIVA_L1O_ENGINE"] = spec["engine"]
    if spec.get("folds"):
        # epang_l1o.py caps K at the number of leaves, so a huge K means "one per fold".
        env["SATIVA_EPANG_FOLDS"] = str(spec["folds"])
    env.update(spec.get("env", {}))
    # Keep thread count honest: only the -T we ask for, no BLAS/OpenMP surprises.
    for var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        env[var] = str(spec["threads"])

    print(f"  [run ] {dataset['name']} / {condition} ...", flush=True)
    meter = [cfg["py3"], str(HERE / "_timeit.py"),
             "--cwd", str(run_dir),
             "--stdout", str(run_dir / "stdout.log")]
    if cfg.get("cpus"):
        command = ["taskset", "-c", cfg["cpus"]] + command
    if timeout:
        meter += ["--timeout", str(timeout)]
    meter += ["--"] + command

    started = datetime.now().isoformat(timespec="seconds")
    proc = subprocess.run(meter, capture_output=True, text=True, env=env)
    if proc.returncode != 0:
        sys.exit(f"timing helper failed:\n{proc.stdout}\n{proc.stderr}")
    measured = json.loads(proc.stdout.strip().splitlines()[-1])

    log_path = run_dir / f"{run_name}.log"
    record = {
        "dataset": dataset["name"],
        "n_seqs": dataset["n_seqs"],
        "aln_len": dataset["aln_len"],
        "condition": condition,
        "rep": rep,
        "impl": spec["impl"],
        "engine": spec["engine"],
        "threads": spec["threads"],
        "folds": spec.get("folds"),
        "env_overrides": spec.get("env", {}),
        "cfg_overrides": spec.get("cfg", {}),
        "started_at": started,
        "host": socket.gethostname(),
        "command": " ".join(command),
        "timeout_sec": timeout,
        "scratch": str(cfg["scratch"]),
        "timing_trusted": cfg["timing_trusted"],
        "cpus": cfg.get("cpus"),
        "n_mislabels": count_mislabels(run_dir / f"{run_name}.mis"),
    }
    record.update(measured)
    record.update(parse_sativa_log(log_path))
    record["status"] = ("timeout" if measured["timed_out"]
                        else "ok" if measured["returncode"] == 0 else "error")

    # Keep only what the analysis needs; the .jplace files reach hundreds of MB.
    out_dir.mkdir(parents=True, exist_ok=True)
    for pattern in ("*.mis", "*.log", "stdout.log"):
        for artefact in run_dir.glob(pattern):
            shutil.copy(artefact, out_dir / artefact.name)
    result_path.write_text(json.dumps(record, indent=2))
    phase = ""
    if "l1out_sec" in record:
        phase = f" (reftree {record['reftree_sec']:.0f}s, l1o {record['l1out_sec']:.0f}s)"
    print(f"        {record['status']}: {record['wall_sec']:.1f}s{phase}, "
          f"peak {record['max_rss_mb']:.0f} MB, mislabels={record['n_mislabels']}", flush=True)

    shutil.rmtree(run_dir, ignore_errors=True)
    shutil.rmtree(temp_dir, ignore_errors=True)

    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--datasets-json", default=str(BENCH_ROOT / "data" / "datasets.json"))
    parser.add_argument("--datasets", default="all", help="comma-separated dataset names, or 'all'")
    parser.add_argument("--conditions", default="all",
                        help="comma-separated condition names, or 'all'. Available: "
                             + ", ".join(CONDITIONS))
    parser.add_argument("--timeout", type=float, default=7200,
                        help="per-run wall-clock cap in seconds (0 = none)")
    parser.add_argument("--runs-dir", default=str(BENCH_ROOT / "results" / "runs"))
    parser.add_argument("--scratch",
                        default=os.environ.get("BENCH_SCRATCH", "/tmp/sativa_speedup_bench"),
                        help="local-disk working root; the project filesystem is too slow to time on")
    parser.add_argument("--force", action="store_true", help="re-run even if run.json exists")
    parser.add_argument("--cpus", default=None,
                        help="pin this batch to a taskset core list (e.g. \"0-7\"). Several "
                             "batches on disjoint core sets run at once without fighting for "
                             "cores; they still share memory bandwidth and last-level cache, "
                             "so pair it with --untimed unless a few percent do not matter")
    parser.add_argument("--untimed", action="store_true",
                        help="this batch answers an agreement question, not a timing one: it "
                             "may share the machine with other batches, and its wall_sec is "
                             "recorded with timing_trusted=false so no figure uses it")
    parser.add_argument("--repeat", type=int, default=1,
                        help="replicates per (dataset, condition); the plot takes the median")
    parser.add_argument("--skip-after-timeout", action="store_true", default=True,
                        help="once a condition times out at size n, skip it for larger n")
    args = parser.parse_args()

    cfg = {
        "py2": env_path("PY2"),
        "py3": env_path("PY3"),
        "py3_env": env_path("PY3_ENV"),
        "sativa_port": env_path("SATIVA_PORT"),
        "sativa_upstream_py2": env_path("SATIVA_UPSTREAM_PY2"),
        "sativa_upstream_py3": env_path("SATIVA_UPSTREAM_PY3"),
        "mode": os.environ.get("SATIVA_MODE", "ultrafast"),
        "conf": os.environ.get("SATIVA_CONF", "0.4"),
        "taxcode": os.environ.get("SATIVA_TAXCODE", "BOT"),
        "seed": os.environ.get("SATIVA_SEED", "42"),
        "runs_dir": args.runs_dir,
        "scratch": args.scratch,
        "timing_trusted": not args.untimed,
        "cpus": args.cpus,
    }
    Path(cfg["scratch"]).mkdir(parents=True, exist_ok=True)

    catalog = json.loads(Path(args.datasets_json).read_text())["datasets"]
    wanted = ([d["name"] for d in catalog] if args.datasets == "all"
              else [s.strip() for s in args.datasets.split(",")])
    datasets = sorted((d for d in catalog if d["name"] in wanted), key=lambda d: d["n_seqs"])

    conditions = (list(CONDITIONS) if args.conditions == "all"
                  else [c.strip() for c in args.conditions.split(",")])
    for condition in conditions:
        if condition not in CONDITIONS:
            sys.exit(f"unknown condition: {condition}")

    Path(cfg["runs_dir"]).mkdir(parents=True, exist_ok=True)
    timeout = args.timeout or None
    gave_up = set()
    started = time.time()

    for dataset in datasets:
        print(f"[{dataset['name']}] {dataset['n_seqs']} seqs x {dataset['aln_len']} cols")
        for condition in conditions:
            if condition in gave_up:
                print(f"  [skip] {condition}: timed out at a smaller size")
                continue
            for rep in range(1, args.repeat + 1):
                record = run_one(cfg, dataset, condition, timeout, force=args.force, rep=rep)
                if record.get("status") == "timeout":
                    if args.skip_after_timeout:
                        gave_up.add(condition)
                    break  # no point repeating a run that ran out of time

    print(f"done in {time.time() - started:.0f}s -- results in {cfg['runs_dir']}")


if __name__ == "__main__":
    main()
