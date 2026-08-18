#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Do the modified and the original SATIVA flag the same sequences?

Reads the `.mis` files each run left behind and compares every condition to a
reference condition on the same dataset:

  matched          same sequence flagged by both
  only_reference   flagged by the reference only  (miss)
  only_test        flagged by the test only       (extra)
  same_rank        of the matched ones, the same MislabeledLevel
  same_proposal    of the matched ones, the same ProposedLabel too

Replicate 1 of each condition is used (SATIVA is deterministic at a fixed seed;
the replicates exist for timing, not for the calls).
"""

import argparse
import json
from pathlib import Path

BENCH_ROOT = Path(__file__).resolve().parent.parent
REFERENCE = "upstream_py3_raxml_T1"


def parse_mis(path):
    """SeqID -> (rank, original, proposed, confidence, orig_path, proposed_path).

    Columns 0-2 are what `compare` uses; the rest is only carried so `--details` can
    name the disagreeing sequences without re-reading the files.
    """
    calls = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(";"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            fields += [""] * (8 - len(fields))
            calls[fields[0]] = (fields[1], fields[2], fields[3], fields[4],
                                fields[5], fields[6])
    return calls


def collect(runs_dir):
    """(dataset, condition) -> calls, for replicate 1 of successful runs."""
    calls = {}
    for run_json in sorted(Path(runs_dir).glob("*/run.json")):
        record = json.loads(run_json.read_text())
        if record["status"] != "ok" or record.get("rep", 1) != 1:
            continue
        mis = list(run_json.parent.glob("*.mis"))
        if not mis:
            continue
        calls[(record["dataset"], record["condition"])] = {
            "n_seqs": record["n_seqs"],
            "calls": parse_mis(mis[0]),
        }
    return calls


def compare(reference, test):
    ref_ids, test_ids = set(reference), set(test)
    matched = ref_ids & test_ids
    same_rank = sum(1 for i in matched if reference[i][0] == test[i][0])
    same_proposal = sum(1 for i in matched
                        if reference[i][0] == test[i][0] and reference[i][2] == test[i][2])
    return {
        "n_reference": len(ref_ids),
        "n_test": len(test_ids),
        "matched": len(matched),
        "only_reference": len(ref_ids - test_ids),
        "only_test": len(test_ids - ref_ids),
        "same_rank": same_rank,
        "same_proposal": same_proposal,
        "recall": len(matched) / len(ref_ids) if ref_ids else None,
        "precision": len(matched) / len(test_ids) if test_ids else None,
        "identical": ref_ids == test_ids and same_proposal == len(matched),
    }


def disagreements(reference, test):
    """One row per sequence the two versions do not call identically."""
    ref_ids, test_ids = set(reference), set(test)
    rows = []
    for seq_id in sorted(ref_ids & test_ids):
        ref, tst = reference[seq_id], test[seq_id]
        if ref[0] != tst[0]:
            kind = "rank differs"
        elif ref[2] != tst[2]:
            kind = "same rank, other proposal"
        else:
            continue
        rows.append({"seq_id": seq_id, "kind": kind,
                     "reference_rank": ref[0], "reference_original": ref[1],
                     "reference_proposed": ref[2], "reference_conf": ref[3],
                     "test_rank": tst[0], "test_original": tst[1],
                     "test_proposed": tst[2], "test_conf": tst[3],
                     "taxonomy": ref[4]})
    for seq_id in sorted(ref_ids - test_ids):
        ref = reference[seq_id]
        rows.append({"seq_id": seq_id, "kind": "flagged by reference only",
                     "reference_rank": ref[0], "reference_original": ref[1],
                     "reference_proposed": ref[2], "reference_conf": ref[3],
                     "test_rank": "", "test_original": "", "test_proposed": "",
                     "test_conf": "", "taxonomy": ref[4]})
    for seq_id in sorted(test_ids - ref_ids):
        tst = test[seq_id]
        rows.append({"seq_id": seq_id, "kind": "flagged by test only",
                     "reference_rank": "", "reference_original": "",
                     "reference_proposed": "", "reference_conf": "",
                     "test_rank": tst[0], "test_original": tst[1],
                     "test_proposed": tst[2], "test_conf": tst[3],
                     "taxonomy": tst[4]})
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", default=str(BENCH_ROOT / "results" / "runs"))
    parser.add_argument("--reference", default=REFERENCE)
    parser.add_argument("--out", default=str(BENCH_ROOT / "results" / "concordance.tsv"))
    parser.add_argument("--details", nargs="?", const=str(BENCH_ROOT / "results" / "concordance_details.tsv"),
                        default=None,
                        help="also write one row per disagreeing sequence (which sequence, "
                             "which rank, which proposal, with what confidence)")
    args = parser.parse_args()

    calls = collect(args.runs_dir)
    datasets = sorted({(key[0], value["n_seqs"]) for key, value in calls.items()},
                      key=lambda item: item[1])

    columns = ["dataset", "n_seqs", "condition", "n_reference", "n_test", "matched",
               "only_reference", "only_test", "same_rank", "same_proposal",
               "recall", "precision", "identical"]
    rows = []
    for dataset, n_seqs in datasets:
        reference = calls.get((dataset, args.reference))
        if reference is None:
            print(f"[warn] no reference run ({args.reference}) for {dataset}, skipped")
            continue
        for (other_dataset, condition), value in sorted(calls.items()):
            if other_dataset != dataset or condition == args.reference:
                continue
            row = {"dataset": dataset, "n_seqs": n_seqs, "condition": condition}
            row.update(compare(reference["calls"], value["calls"]))
            rows.append(row)

    with open(args.out, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(
                "" if row.get(c) is None
                else f"{row[c]:.3f}" if isinstance(row.get(c), float)
                else str(row.get(c, "")) for c in columns) + "\n")

    for row in rows:
        flag = "identical" if row["identical"] else (
            f"recall {row['recall']:.2f} / precision {row['precision']:.2f}"
            if row["recall"] is not None else "n/a")
        print(f"{row['dataset']:>7} {row['condition']:<24} "
              f"{row['n_reference']:>5} vs {row['n_test']:>5} flagged -> {flag}")
    print(f"\nwrote {args.out} ({len(rows)} comparisons, reference = {args.reference})")

    if args.details:
        detail_columns = ["dataset", "n_seqs", "condition", "seq_id", "kind",
                          "reference_rank", "reference_original", "reference_proposed",
                          "reference_conf", "test_rank", "test_original", "test_proposed",
                          "test_conf", "taxonomy"]
        detail_rows = []
        for dataset, n_seqs in datasets:
            reference = calls.get((dataset, args.reference))
            if reference is None:
                continue
            for (other_dataset, condition), value in sorted(calls.items()):
                if other_dataset != dataset or condition == args.reference:
                    continue
                for row in disagreements(reference["calls"], value["calls"]):
                    row.update({"dataset": dataset, "n_seqs": n_seqs, "condition": condition})
                    detail_rows.append(row)
        with open(args.details, "w") as handle:
            handle.write("\t".join(detail_columns) + "\n")
            for row in detail_rows:
                handle.write("\t".join(str(row.get(c, "")) for c in detail_columns) + "\n")
        print(f"wrote {args.details} ({len(detail_rows)} disagreeing sequences)")


if __name__ == "__main__":
    main()
