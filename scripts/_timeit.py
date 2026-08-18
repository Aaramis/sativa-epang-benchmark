#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Run one command, report wall time and peak child RSS as JSON on stdout.

A fresh process per measured run, so `RUSAGE_CHILDREN.ru_maxrss` covers this
run only (there is no /usr/bin/time on this host). The command runs in its own
process group, so a timeout kills the whole RAxML/EPA-ng subtree, not just the
python parent.
"""

import argparse
import json
import os
import resource
import signal
import subprocess
import sys
import time


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--timeout", type=float, default=None)
    parser.add_argument("--cwd", default=None)
    parser.add_argument("--stdout", default=None)
    parser.add_argument("cmd", nargs=argparse.REMAINDER)
    args = parser.parse_args()

    cmd = args.cmd[1:] if args.cmd and args.cmd[0] == "--" else args.cmd
    if not cmd:
        sys.exit("no command given")

    out = open(args.stdout, "wb") if args.stdout else subprocess.DEVNULL
    start = time.time()
    timed_out = False
    proc = subprocess.Popen(cmd, cwd=args.cwd, stdout=out, stderr=subprocess.STDOUT,
                            start_new_session=True)
    try:
        rc = proc.wait(timeout=args.timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        os.killpg(proc.pid, signal.SIGKILL)
        rc = proc.wait()
    wall = time.time() - start
    if out is not subprocess.DEVNULL:
        out.close()

    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    print(json.dumps({
        "wall_sec": round(wall, 2),
        "returncode": rc,
        "timed_out": timed_out,
        "max_rss_mb": round(usage.ru_maxrss / 1024.0, 1),   # ru_maxrss is KB on Linux
        "user_cpu_sec": round(usage.ru_utime, 2),
        "sys_cpu_sec": round(usage.ru_stime, 2),
    }))


if __name__ == "__main__":
    main()
