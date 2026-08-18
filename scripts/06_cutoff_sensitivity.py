#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Does raising SATIVA's confidence cutoff make the two engines agree?

Every disagreement we could name -- the Ilex, the Bouteloua pair, the divergent phylum
proposals -- sits between 0.4 and 0.6 confidence. Giorgos' COX1 workflow keeps only calls
at genus level or above with confidence >= 0.9. This re-reads the .mis files already on
disk, applies a stricter cutoff to both sides, and recomputes the agreement, so the
question is answered without running anything again.

The runs were made with -C 0.4, so cutoffs below that cannot be explored here; the -C 0
probe (scripts/18) is what opens the other direction.
"""

import argparse
import glob
import json
from pathlib import Path

BENCH_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = "upstream_py3_raxml_T1"
RANK_ORDER = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def parse_mis(path, min_conf=0.0, max_rank=None):
    """SeqID -> (rank, proposed) for calls passing a confidence / rank filter."""
    calls = {}
    for line in open(path):
        if not line.strip() or line.startswith(";"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        try:
            conf = float(f[4])
        except ValueError:
            continue
        if conf < min_conf:
            continue
        if max_rank is not None and f[1] in RANK_ORDER \
                and RANK_ORDER.index(f[1]) > RANK_ORDER.index(max_rank):
            continue
        calls[f[0]] = (f[1], f[3])
    return calls


def collect(runs_dir, min_conf, max_rank):
    out = {}
    for run_json in sorted(Path(runs_dir).glob("*/run.json")):
        rec = json.loads(run_json.read_text())
        if rec["status"] != "ok" or rec.get("rep", 1) != 1:
            continue
        mis = glob.glob(str(run_json.parent / "*.mis"))
        if mis:
            out[(rec["dataset"], rec["condition"])] = (rec["n_seqs"],
                                                       parse_mis(mis[0], min_conf, max_rank))
    return out


def by_confidence_bin(runs_dir, reference, conditions, bins):
    """Does the reference's own confidence predict whether a call is reproduced?

    A cutoff is only worth raising if low confidence calls are the unstable ones. For each
    band of the reference's confidence, this counts how many of its calls the other version
    reproduces at the same rank, pooled over every alignment size.
    """
    import collections
    calls = {}
    for run_json in sorted(Path(runs_dir).glob("*/run.json")):
        rec = json.loads(run_json.read_text())
        if rec["status"] != "ok" or rec.get("rep", 1) != 1:
            continue
        mis = glob.glob(str(run_json.parent / "*.mis"))
        if mis:
            calls[(rec["dataset"], rec["condition"])] = parse_mis_with_conf(mis[0])

    for condition in conditions:
        tally = collections.defaultdict(lambda: [0, 0])
        for (dataset, cond), reference_calls in calls.items():
            if cond != reference:
                continue
            test = calls.get((dataset, condition))
            if not test:
                continue
            for seq, (rank, conf) in reference_calls.items():
                for lo, hi in bins:
                    if lo <= conf < hi:
                        tally[(lo, hi)][1] += 1
                        if seq in test and test[seq][0] == rank:
                            tally[(lo, hi)][0] += 1
        print(f"\n{condition}: reference calls reproduced, by the reference's own confidence")
        for band in bins:
            got, total = tally[band]
            if total:
                print(f"   {band[0]:.2f} to {min(band[1], 1.0):.2f}: "
                      f"{got:>4}/{total:<4} = {got / total:.2f}")


def parse_mis_with_conf(path):
    """SeqID -> (rank, confidence)."""
    calls = {}
    for line in open(path):
        if line.startswith(";") or not line.strip():
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 5:
            continue
        try:
            calls[f[0]] = (f[1], float(f[4]))
        except ValueError:
            continue
    return calls


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--runs-dir", default=str(BENCH_ROOT / "results" / "runs"))
    parser.add_argument("--reference", default=REFERENCE)
    parser.add_argument("--conditions", default="port_py3_epang_T8,port_py3_epang_T8_K25")
    parser.add_argument("--by-bin", action="store_true",
                        help="instead of cutoffs, report reproduction rate per confidence band")
    parser.add_argument("--cutoffs", default="0.4,0.6,0.7,0.8,0.9,0.95")
    parser.add_argument("--max-rank", default=None,
                        help="also drop calls below this rank, e.g. Genus")
    args = parser.parse_args()

    conditions = [c.strip() for c in args.conditions.split(",")]
    print(f"reference = {args.reference}"
          + (f", ranks kept: {args.max_rank} and above" if args.max_rank else ""))

    if args.by_bin:
        by_confidence_bin(args.runs_dir, args.reference, conditions,
                          [(0.4, 0.5), (0.5, 0.7), (0.7, 0.9), (0.9, 1.0001)])
        return

    header = f"{'dataset':>7} {'condition':<26} {'cutoff':>6} {'ref':>5} {'ours':>5} " \
             f"{'common':>6} {'recall':>7} {'prec':>6}"
    print(header)
    for cutoff in [float(c) for c in args.cutoffs.split(",")]:
        calls = collect(args.runs_dir, cutoff, args.max_rank)
        datasets = sorted({(ds, n) for (ds, _), (n, _) in calls.items()}, key=lambda x: x[1])
        for dataset, _ in datasets:
            ref = calls.get((dataset, args.reference))
            if not ref:
                continue
            for condition in conditions:
                ours = calls.get((dataset, condition))
                if not ours:
                    continue
                ref_ids, our_ids = set(ref[1]), set(ours[1])
                common = ref_ids & our_ids
                recall = len(common) / len(ref_ids) if ref_ids else float("nan")
                prec = len(common) / len(our_ids) if our_ids else float("nan")
                print(f"{dataset:>7} {condition:<26} {cutoff:>6.2f} {len(ref_ids):>5} "
                      f"{len(our_ids):>5} {len(common):>6} {recall:>7.3f} {prec:>6.3f}")
        print()


if __name__ == "__main__":
    main()
