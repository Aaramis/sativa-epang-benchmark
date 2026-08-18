#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Ground-truth test of the failure mode Daniel Lundin describes.

His observation on GTDB: a handful of S. aureus sequences sit inside E. coli and SATIVA
does not flag them. His explanation -- the leave-one-out removes one sequence at a time,
so the four remaining copies of the same wrong label anchor the placement and confirm the
error. If that is right, a K-fold leave-one-out should do *better* than the strict one,
because a fold takes several sequences out at once.

This builds the experiment instead of arguing about it: take a real subset, relabel k
sequences of one species with the name of a distant genus, and see which method flags
them. k = 1 is the case the strict leave-one-out handles; k = 2, 3, 5 are the cases it is
supposed to be blind to.

Usage:
    07_injected_mislabels.py build --k 1,2,3,5      # writes data_injected/
    07_injected_mislabels.py score                   # reads the runs back
"""

import argparse
import json
import random
from pathlib import Path

BENCH_ROOT = Path(__file__).resolve().parent.parent
RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
RANK_TITLES = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]


def read_tax(path):
    out = {}
    for line in open(path):
        sid, tax = line.rstrip("\n").split("\t")
        out[sid] = tax.split(";")
    return out


def build(args):
    source = Path(args.source)
    tax = read_tax(source / f"{source.name}.tax")
    fasta = (source / f"{source.name}.fasta").read_text()

    by_species = {}
    for sid, ranks in tax.items():
        by_species.setdefault(";".join(ranks), []).append(sid)

    # Donor: a species with enough members to relabel; recipient: a different family, so
    # the injected error is visible at genus level and above.
    donors = sorted((s for s, ids in by_species.items() if len(ids) >= max(args.k_values)),
                    key=lambda s: (-len(by_species[s]), s))
    if not donors:
        raise SystemExit("no species with enough members in this subset")
    if args.donor_index >= len(donors):
        raise SystemExit(f"only {len(donors)} species large enough; donor-index too high")
    donor = donors[args.donor_index]
    donor_ranks = donor.split(";")
    recipients = [s.split(";") for s in by_species
                  if s.split(";")[4] != donor_ranks[4] and s.split(";")[5] != "NA"]
    recipient = sorted(recipients, key=lambda r: ";".join(r))[args.donor_index % len(recipients)]

    print(f"donor species : {donor}  ({len(by_species[donor])} sequences)")
    print(f"wrong label   : {';'.join(recipient)}")

    rng = random.Random(args.seed)
    members = sorted(by_species[donor])
    rng.shuffle(members)

    outroot = Path(args.outdir)
    manifest = []
    for k in args.k_values:
        injected = sorted(members[:k])
        name = f"{source.name}_d{args.donor_index}inj{k}" if args.donor_index else f"{source.name}_inj{k}"
        dest = outroot / name
        dest.mkdir(parents=True, exist_ok=True)
        (dest / f"{name}.fasta").write_text(fasta)
        with open(dest / f"{name}.tax", "w") as handle:
            for sid, ranks in tax.items():
                out = list(recipient) if sid in injected else ranks
                handle.write(f"{sid}\t{';'.join(out)}\n")
        (dest / "injected.json").write_text(json.dumps(
            {"injected": injected, "donor": donor, "wrong_label": ";".join(recipient)},
            indent=2))
        manifest.append({"name": name, "n_seqs": len(tax), "aln_len": 0,
                         "fasta": str(dest / f"{name}.fasta"),
                         "tax": str(dest / f"{name}.tax"), "k_injected": k})
        print(f"  {name}: {k} sequence(s) relabelled")

    catalog = outroot / "datasets.json"
    existing = json.loads(catalog.read_text())["datasets"] if catalog.exists() else []
    keep = [e for e in existing if e["name"] not in {m["name"] for m in manifest}]
    catalog.write_text(json.dumps(
        {"source": str(source), "seed": args.seed, "datasets": keep + manifest}, indent=2))
    print(f"wrote {outroot / 'datasets.json'}")


def score(args):
    catalog = json.loads((Path(args.outdir) / "datasets.json").read_text())["datasets"]
    runs = Path(args.runs_dir)
    print(f"{'dataset':<14} {'k':>2} {'condition':<26} {'found':>6} {'at rank':<10} {'total flags':>11}")
    for entry in catalog:
        info = json.loads((Path(entry["tax"]).parent / "injected.json").read_text())
        injected = set(info["injected"])
        for run_json in sorted(runs.glob(f"{entry['name']}__*/run.json")):
            rec = json.loads(run_json.read_text())
            if rec["status"] != "ok":
                continue
            mis = list(run_json.parent.glob("*.mis"))
            if not mis:
                continue
            found, ranks, total = set(), [], 0
            for line in open(mis[0]):
                if not line.strip() or line.startswith(";"):
                    continue
                f = line.rstrip("\n").split("\t")
                total += 1
                if f[0] in injected:
                    found.add(f[0])
                    ranks.append(f[1])
            print(f"{entry['name']:<14} {entry['k_injected']:>2} {rec['condition']:<26} "
                  f"{len(found):>3}/{len(injected):<2} {','.join(sorted(set(ranks))) or '-':<10} "
                  f"{total:>11}")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)

    b = sub.add_parser("build")
    b.add_argument("--source", default=str(BENCH_ROOT / "data_ascii" / "n400"))
    b.add_argument("--outdir", default=str(BENCH_ROOT / "data_injected"))
    b.add_argument("--k", dest="k_values", default="1,2,3,5",
                   type=lambda s: [int(x) for x in s.split(",")])
    b.add_argument("--seed", type=int, default=42)
    b.add_argument("--donor-index", type=int, default=0,
                   help="0 = the largest species, 1 = the next one, ... (a different "
                        "donor/recipient pair, so the finding does not rest on one clade)")
    b.set_defaults(func=build)

    s = sub.add_parser("score")
    s.add_argument("--outdir", default=str(BENCH_ROOT / "data_injected"))
    s.add_argument("--runs-dir", default=str(BENCH_ROOT / "results" / "runs_injected"))
    s.set_defaults(func=score)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
