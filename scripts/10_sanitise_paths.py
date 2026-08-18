#!/usr/bin/env python3
# Part of the SATIVA EPA-ng benchmark. What the numbers mean: ../RESULTS.md
"""Replace machine specific paths in the recorded runs with placeholders.

Every run.json, .log and .mis under results/ carries the command line that produced it,
which spells out the interpreter, the working directory and the scratch directory of the
machine it ran on. That is noise to a reader and it names directories that have nothing to
do with this work, so it is rewritten before publication:

    <conda env>/bin/python            -> $PY3
    <bench root>/                     -> ./
    /tmp/<scratch>/                   -> $SCRATCH/

Idempotent: running it twice changes nothing. It rewrites files in place, so run it on a
copy if the originals matter to you.
"""

import argparse
import re
from pathlib import Path

BENCH_ROOT = Path(__file__).resolve().parent.parent


def substitutions(bench_root, py3_env, scratch):
    rules = []
    if py3_env:
        rules.append((re.escape(str(py3_env).rstrip("/")) + r"/bin/python\b", "$PY3"))
        rules.append((re.escape(str(py3_env).rstrip("/")), "$PY3_ENV"))
    rules.append((re.escape(str(bench_root).rstrip("/")) + "/", "./"))
    if scratch:
        rules.append((re.escape(str(scratch).rstrip("/")), "$SCRATCH"))
    # Anything left that points into a home or a mounted project tree.
    rules.append((r"/pools/[^\s\"']*/(?=[a-zA-Z0-9_.-]+)", ""))
    rules.append((r"/home/[a-z0-9_-]+/[^\s\"']*/(?=[a-zA-Z0-9_.-]+)", ""))
    return [(re.compile(pattern), replacement) for pattern, replacement in rules]


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default=str(BENCH_ROOT / "results"))
    parser.add_argument("--bench-root", default=str(BENCH_ROOT))
    parser.add_argument("--py3-env", default=None,
                        help="conda env whose python ran the benchmark")
    parser.add_argument("--scratch", default="/tmp/sativa_speedup_bench")
    parser.add_argument("--also", nargs="*", default=[],
                        help="extra files outside --root to rewrite")
    args = parser.parse_args()

    rules = substitutions(args.bench_root, args.py3_env, args.scratch)
    targets = [p for p in Path(args.root).rglob("*")
               if p.is_file() and p.suffix in {".json", ".log", ".mis", ".tsv"}]
    targets += [Path(p) for p in args.also]

    changed = 0
    for path in targets:
        try:
            text = path.read_text()
        except UnicodeDecodeError:
            continue
        new = text
        for pattern, replacement in rules:
            new = pattern.sub(replacement, new)
        if new != text:
            path.write_text(new)
            changed += 1
    print(f"{changed} of {len(targets)} files rewritten")


if __name__ == "__main__":
    main()
