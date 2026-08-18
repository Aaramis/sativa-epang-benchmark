#!/usr/bin/env python3
"""Aggregate results/runs/*/run.json into summary.tsv and draw the figure.

One figure, two questions:

  A  where the time goes -- reference tree (shared, RAxML) vs leave-one-out
     (the step we replaced), one panel per alignment size
  B  do we flag the same sequences as unmodified SATIVA?

Replicates are collapsed with the median. A run that hit the wall-clock cap is drawn as such -- "no result at the
cap" is the result, not missing data.
"""

import argparse
import importlib.util
import json
import textwrap
from pathlib import Path
from statistics import median

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch

BENCH_ROOT = Path(__file__).resolve().parent.parent

# Categorical slots of the validated palette (light mode) + neutral ink.
BLUE, ORANGE, AQUA, YELLOW = "#2a78d6", "#eb6834", "#1baf7a", "#eda100"
INK, INK_SOFT, MUTED, SHARED, GRID = "#0b0b0b", "#52514e", "#8a8983", "#c9c8c2", "#ececea"
# Everything that is neither the reference tree nor the leave-one-out: SATIVA's own
# bookkeeping (refjson, the final EPA check, I/O) plus whatever the whole-second phase
# timer rounded away. It is not part of either engine, so it wears neither engine colour.
RESIDUAL = "#e4e3de"

BASELINE = "upstream_py3_raxml_T1"          # unmodified SATIVA v0.9.3, python 3
# Panel B rows. Thread count does not change the calls (T1, T2, T4, T8 and T16 give
# byte-identical .mis files), so the number of folds is the only axis worth showing.
AGREEMENT_ROWS = [
    ("port_py3_epang_T8",     "EPA-ng, K=5 folds"),
    ("port_py3_epang_T8_K25", "EPA-ng, K=25 folds"),
]
CONTROL_ROW = "upstream_py2_raxml_T1"       # unmodified SATIVA, python-2 build

TITLE = "SATIVA with EPA-ng: runtime and agreement"

# "Original" is upstream v0.9.3 itself -- the code ours is a three-hunk diff against --
# so the labels name it rather than saying "original", which invites the wrong baseline.
TIME_BARS = [
    ("upstream_py3_raxml_T1", "upstream v0.9.3 · 1 thread",       ORANGE),
    ("upstream_py3_raxml_T8", "upstream v0.9.3 · 8 threads",      ORANGE),
    ("port_py3_epang_T1",     "EPA-ng, K=5 folds · 1 thread",     AQUA),
    ("port_py3_epang_T8",     "EPA-ng, K=5 folds · 8 threads",    AQUA),
    ("port_py3_epang_T1_K25", "EPA-ng, K=25 folds · 1 thread",    AQUA),
    ("port_py3_epang_T8_K25", "EPA-ng, K=25 folds · 8 threads",   AQUA),
]

# SATIVA logs its phase split in whole seconds. Below a few hundred sequences the
# reference tree takes under a second and the split degenerates to "0 s reftree",
# so panel A starts where the numbers mean something. Every size stays in
# summary.tsv and in panel B, where the resolution does not matter.
MIN_SIZE_TIME_PANEL = 400

# SATIVA changes what it asks RAxML to do as the tree grows (epac/config.py):
# above CAT_GAMMA_THRES=500 taxa the model drops from GTRGAMMA to GTRCAT, and above
# EPA_HEUR_THRES=1000 it turns on the EPA heuristic, so only a fraction of the branches
# get a thorough insertion (-G 0.5*1000/n). EPA-ng keeps the full GTRGAMMA computation on
# every branch at every size. Past 500 taxa the orange bar is therefore a cheaper
# approximation, not the same computation shortened -- which is exactly why the ratio
# stops growing there. The badge says so on the panel where it happens.
def raxml_regime(n_seqs):
    if n_seqs > 1000:
        return f"RAxML here: GTRCAT + EPA heuristic (-G, ~{50000 / n_seqs:.0f} % of branches)"
    if n_seqs > 500:
        return "RAxML here: GTRCAT, every branch"
    return "RAxML here: GTRGAMMA, every branch"


# --------------------------------------------------------------------------- data

def load_runs(runs_dir):
    return [json.loads(p.read_text()) for p in sorted(Path(runs_dir).glob("*/run.json"))]


def aggregate(records):
    """(dataset, condition) -> median over replicates, with the min-max spread."""
    buckets = {}
    for rec in records:
        buckets.setdefault((rec["dataset"], rec["condition"]), []).append(rec)
    table = {}
    for key, reps in buckets.items():
        ok = [r for r in reps if r["status"] == "ok"]
        row = dict(reps[0])
        row["n_reps"] = len(ok) or len(reps)
        row["status"] = "ok" if ok else reps[0]["status"]
        if ok:
            for field in ("wall_sec", "reftree_sec", "l1out_sec", "max_rss_mb"):
                values = [r[field] for r in ok if r.get(field) is not None]
                row[field] = median(values) if values else None
            walls = [r["wall_sec"] for r in ok]
            row["wall_min"], row["wall_max"] = min(walls), max(walls)
        table[key] = row
    return table


def write_summary(table, out_path):
    columns = ["dataset", "n_seqs", "aln_len", "condition", "impl", "engine", "threads",
               "status", "n_reps", "wall_sec", "wall_min", "wall_max", "reftree_sec",
               "l1out_sec", "max_rss_mb", "n_mislabels"]
    rows = sorted(table.values(), key=lambda r: (r["n_seqs"], r["condition"]))
    with open(out_path, "w") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join("" if row.get(c) is None else str(row.get(c, ""))
                                   for c in columns) + "\n")
    return rows


def load_concordance(runs_dir, reference, min_conf=0.0):
    """Reuse 05_concordance.py rather than re-implementing the .mis comparison.

    `min_conf` filters both sides at the same confidence before comparing, which is how
    panel C is drawn: the disagreements the pipeline would actually integrate depend on
    where the cutoff sits, not only on the engine.
    """
    spec = importlib.util.spec_from_file_location(
        "concordance", Path(__file__).resolve().parent / "05_concordance.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    calls = module.collect(runs_dir)
    if min_conf > 0:
        for value in calls.values():
            value["calls"] = {seq: call for seq, call in value["calls"].items()
                              if call[3] and float(call[3]) >= min_conf}
    rows = []
    for (dataset, condition), value in sorted(calls.items()):
        if condition == reference:
            continue
        ref = calls.get((dataset, reference))
        if ref is None:
            continue
        row = {"dataset": dataset, "n_seqs": value["n_seqs"], "condition": condition}
        row.update(module.compare(ref["calls"], value["calls"]))
        rows.append(row)
    return rows


def datasets_in(table):
    seen = {}
    for row in table.values():
        seen[row["dataset"]] = row["n_seqs"]
    return sorted(seen.items(), key=lambda kv: kv[1])


# ------------------------------------------------------------------------ drawing

def fmt_time(seconds):
    if seconds < 90:
        return f"{seconds:.0f} s"
    if seconds < 5400:
        return f"{seconds / 60:.0f} min"
    return f"{seconds / 3600:.1f} h"


def fmt_n(value):
    return f"{value:,}".replace(",", " ")


def bare_axes(ax, xgrid=True):
    ax.set_facecolor("white")
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color(GRID)
    ax.tick_params(colors=INK_SOFT, labelsize=8, length=0)
    if xgrid:
        ax.grid(axis="x", color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def draw_time_panel(ax, table, dataset, n_seqs, timeout_hint, show_labels):
    # Every condition always keeps its row, so a version that has not run yet reads
    # as a gap instead of silently disappearing from the comparison.
    rows = [table.get((dataset, condition)) for condition, _, _ in TIME_BARS]
    labels = [label for _, label, _ in TIME_BARS]
    colors = [color for _, _, color in TIME_BARS]
    if not any(row is not None for row in rows):
        ax.set_axis_off()
        return

    measured = [r for r in rows if r is not None]
    cap = max((r.get("timeout_sec") or timeout_hint) for r in measured)
    finished = [r["wall_sec"] for r in measured if r["status"] == "ok"]
    timed_out = any(r["status"] != "ok" for r in measured)
    xmax = max(finished + ([cap] if timed_out else []))

    for y, (row, color) in enumerate(zip(rows, colors)):
        if row is not None and not row.get("timing_trusted", True):
            row = None  # measured while sharing the machine: good for the calls, not for a time
        if row is None:
            ax.annotate("not run", (0, y), xytext=(4, 0), textcoords="offset points",
                        va="center", fontsize=8, color=MUTED, style="italic")
        elif row["status"] == "ok":
            reftree = row.get("reftree_sec") or 0
            l1out = row.get("l1out_sec") or 0
            other = max(row["wall_sec"] - reftree - l1out, 0)
            # Drawn in the order the run spends the time: SATIVA's own bookkeeping,
            # then the reference tree, then the leave-one-out. Both shared (grey)
            # phases sit together on the left, so the coloured part -- the only thing
            # that differs between versions -- is what the eye compares.
            if other > 0.02 * xmax:
                ax.barh(y, other, height=0.62, color=RESIDUAL, zorder=3)
            else:
                other = 0
            ax.barh(y, reftree, left=other, height=0.62, color=SHARED, zorder=3)
            ax.barh(y, l1out, left=other + reftree, height=0.62, color=color, zorder=3)
            # The replicate spread stays in summary.tsv; drawing it here only added
            # clutter to bars whose message is the median.
            ax.annotate(fmt_time(row["wall_sec"]), (row["wall_sec"], y), xytext=(7, 0),
                        textcoords="offset points", va="center", fontsize=8.5, color=INK)
        else:
            ax.barh(y, cap, height=0.62, color="white", edgecolor=MUTED,
                    hatch="////", linewidth=0.9, zorder=3)
            ax.annotate(f"no result at {fmt_time(cap)}", (cap * 0.5, y), ha="center",
                        va="center", fontsize=8, color=INK,
                        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                                  edgecolor="none"))

    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(labels if show_labels else [""] * len(rows),
                       fontsize=8.5, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, xmax * 1.3)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_title(f"{fmt_n(n_seqs)} sequences", fontsize=10, color=INK, loc="left", pad=22)
    ax.annotate(raxml_regime(n_seqs), (0, 1), xycoords="axes fraction",
                xytext=(0, 8), textcoords="offset points", fontsize=7.8,
                color=MUTED, ha="left", va="bottom")
    ax.set_xlabel("seconds", fontsize=8, color=MUTED)
    ax.xaxis.set_major_locator(plt.MaxNLocator(4))
    bare_axes(ax)


def agreement_entries(concordance_rows, min_size=0):
    """Rows of panel B: one per size, our engine against unmodified SATIVA v0.9.3.

    The python-2 comparison (upstream v0.9.3 against upstream's own last python-2 commit)
    is deliberately NOT drawn: it is a second, independent comparison against the same
    reference, and side by side it read as if ours were being compared to python 2. It
    stays in concordance.tsv / concordance_details.tsv, where it can be quoted as a scale
    for how many calls move under a change upstream did not consider behavioural.
    """
    by_key = {(row["n_seqs"], row["condition"]): row for row in concordance_rows}
    # Same size floor as panel A, so the two panels describe the same runs.
    sizes = sorted({n for n, _ in by_key if n >= min_size})

    entries = []
    for group, n_seqs in enumerate(sizes):
        for condition, label in AGREEMENT_ROWS:
            row = by_key.get((n_seqs, condition))
            if row:
                entries.append((f"{fmt_n(n_seqs)} seqs · {label}", row, False, group))
    return entries


def draw_agreement_panel(ax, entries):
    if not entries:
        ax.set_axis_off()
        return

    ticks, xmax, y, previous_group = [], 0, 0.0, None
    for _, row, _, group in entries:
        if previous_group is not None and group != previous_group:
            y += 0.55                       # a breath between two alignment sizes
        previous_group = group
        ticks.append(y)
        segments = ((row["same_rank"], AQUA),
                    (row["matched"] - row["same_rank"], YELLOW),
                    (row["only_reference"], SHARED),
                    (row["only_test"], BLUE))
        left = 0
        for value, color in segments:
            if value:
                ax.barh(y, value, left=left, height=0.6, color=color, zorder=3)
                left += value
        xmax = max(xmax, left)
        ax.annotate(f"{row['same_rank']}/{row['n_reference']} identical",
                    (left, y), xytext=(8, 0), textcoords="offset points",
                    va="center", fontsize=8.5, color=INK)
        y += 1

    ax.set_yticks(ticks)
    ax.set_yticklabels([label for label, _, _, _ in entries], fontsize=8.5, color=INK)
    for tick, (_, _, is_control, _) in zip(ax.get_yticklabels(), entries):
        if is_control:
            tick.set_color(MUTED)
    ax.set_ylim(ticks[-1] + 0.8, ticks[0] - 0.8)
    ax.set_xlim(0, xmax * 1.35)
    ax.set_xlabel("sequences flagged as mislabelled", fontsize=8, color=MUTED)
    ax.xaxis.set_major_locator(plt.MaxNLocator(6))
    bare_axes(ax)


def wrap(text, fig_width_in, fontsize):
    """Wrap to the canvas: ~0.0072 in per character per point of font size."""
    chars = max(int(fig_width_in / (0.0072 * fontsize)), 20)
    return textwrap.wrap(text, width=chars)


def figure_main(table, concordance_rows, strict_rows, outdir, timeout_hint,
                min_size, strict_cutoff):
    """One figure, laid out in inches from the top down.

    matplotlib's automatic layout cannot see free-floating text, so titles, section
    labels and legends collided with the panels. Here every block declares its height,
    the figure height is their sum, and each block is placed at an explicit offset.
    """
    sizes = [(name, n) for name, n in datasets_in(table) if n >= min_size]
    if not sizes:
        # Progress view: while the benchmark is still on the small sizes, draw what exists
        # rather than leaving a stale file from an earlier run on disk.
        sizes = datasets_in(table)
    if not sizes:
        print("[warn] no dataset measured yet; nothing to draw")
        return
    smallest = sizes[0][1]

    ncols = len(sizes)
    width = max(3.6 * ncols + 2.3, 11.0)
    left_in, right_in, col_gap_in = 2.35, 0.35, 0.75

    legend_cols = 4 if width >= 13.5 else 2
    legend_rows = 4 // legend_cols

    agreement_rows = agreement_entries(concordance_rows, smallest)
    strict_agreement_rows = agreement_entries(strict_rows, smallest)
    subtitle = ("Nested subsets of one ITS alignment. Every run applies the same SATIVA decision rule and builds "
                "its reference tree with RAxML 8.2.3 (grey). What changes is the placement engine: the "
                "leave-one-out (coloured) and, inside \"other steps\", the final confirmation pass.")
    footer = ("Medians of the replicates (spread in summary.tsv). Every run works on local disk: on the "
              "project's network filesystem the I/O dominates and adds several-fold noise. Panel A starts at "
              f"{smallest} sequences and panel B follows it. Below about 400, SATIVA's whole-second phase timer "
              "cannot resolve the split, and the smaller sizes stay in summary.tsv and concordance.tsv. "
              "In panel B the k-fold leave-one-out is an approximation of the per-sequence one, so it drops and "
              "adds a few flags.")
    note_a = ("Upstream v0.9.3 = SATIVA as its authors ship it (python 3; the python-3 port is theirs, "
              "commit 259186e). The EPA-ng version is that same code with the placement engine swapped: the "
              "decision-logic files are byte-identical, the diff is two hunks in sativa.py plus one new module. "
              "Mind the badge "
              "over each panel: past 500 taxa SATIVA drops its own RAxML from GTRGAMMA to GTRCAT, and past 1000 "
              "taxa it thoroughly places on a fraction of the branches only, while EPA-ng keeps the full "
              "computation everywhere. From n=800 on, the orange bar is a cheaper approximation rather than the "
              "same work done faster. That is why the ratio stops growing, and why at 5 402 sequences the "
              "speed-up comes from threads: EPA-ng parallelises over queries, RAxML over 242 alignment columns, "
              "where there is nothing to split.")
    note_c = (f"The same comparison with both sides filtered at confidence >= {strict_cutoff:g}. It does not "
              "make the two agree: with K=25, recall is 0.93 at n=800, 0.86 at n=1 600, 0.81 at n=5 402. And "
              "confidence predicts reproducibility only at the bottom of its range: calls between 0.40 and 0.50 "
              "are reproduced 0.69 of the time, against about 0.8 above, with a call at 0.95 no safer than one "
              "at 0.6.")
    note_b = ("Each row: the EPA-ng version against unmodified SATIVA v0.9.3, sequence by sequence, cutoff 0.4 "
              "on both sides. K is the number of folds the leave-one-out is split into, the only approximation "
              "in the modification and the knob that closes most of the gap. At a fixed seed the replicates give "
              "identical calls, so none of this is run-to-run noise.")
    note_a = ("Upstream v0.9.3 is SATIVA as its authors ship it, python 3 included. The EPA-ng version is that "
              "same code with the placement engine swapped: the decision-logic files are byte-identical, the "
              "diff is two hunks in sativa.py, one in config.py and one new module. The badge over each panel "
              "says what SATIVA asks RAxML to do at that size: past 500 taxa GTRCAT instead of GTRGAMMA, past "
              "1 000 taxa a thorough insertion on a fraction of the branches only. Those shortcuts barely change "
              "its output (119 flags against 118 at n=1 600) but they do make the orange bar a cheaper "
              "computation from n=800 on, which is why the ratio stops growing. The EPA-ng version keeps its own "
              "reference tree on GTRGAMMA at every size, so that EPA-ng reads a fitted shape parameter.")
    subtitle_lines = wrap(subtitle, width - 0.2, 9.2)
    note_a_lines = wrap(note_a, width - 0.2, 8.6)
    note_b_lines = wrap(note_b, width - 0.2, 8.6)
    note_c_lines = wrap(note_c, width - 0.2, 8.6)
    footer_lines = wrap(footer, width - 0.2, 8.0)

    # A block is the axes plus the room its own decorations need:
    # PANEL_TITLE_IN above (the "400 sequences" line), AXIS_IN below (ticks + label).
    panel_title_in, axis_in = 0.62, 0.62
    panel_a_in = 0.62 * len(TIME_BARS) + panel_title_in + axis_in
    def panel_height(entries):
        groups = len({entry[3] for entry in entries}) if entries else 1
        return (0.42 * max(len(entries), 1) + 0.23 * max(groups - 1, 0) + 0.12 + axis_in)

    panel_b_in = panel_height(agreement_rows)
    panel_c_in = panel_height(strict_agreement_rows)
    blocks = [
        ("margin", 0.30),
        ("title", 0.34),
        ("subtitle", 0.20 * len(subtitle_lines) + 0.20),
        ("section_a", 0.34),
        ("section_a_note", 0.19 * len(note_a_lines) + 0.13),
        ("panel_a", panel_a_in),
        ("legend_a", 0.30 * legend_rows + 0.06),
        ("section_b", 0.34),
        ("section_b_note", 0.19 * len(note_b_lines) + 0.13),
        ("panel_b", panel_b_in),
        ("legend_b", 0.30 * legend_rows + 0.06),
        ("section_c", 0.34),
        ("section_c_note", 0.19 * len(note_c_lines) + 0.13),
        ("panel_c", panel_c_in),
        ("footer", 0.17 * len(footer_lines) + 0.30),
    ]
    height = sum(h for _, h in blocks)
    fig = plt.figure(figsize=(width, height), facecolor="white")

    def y_of(name):
        """Top edge of a block, as a figure fraction."""
        offset = 0.0
        for block, block_h in blocks:
            if block == name:
                return 1 - offset / height
            offset += block_h
        raise KeyError(name)

    def block_bottom(name):
        block_h = dict(blocks)[name]
        return y_of(name) - block_h / height

    left_frac, x_right = left_in / width, 1 - right_in / width
    col_w = (x_right - left_frac - (ncols - 1) * col_gap_in / width) / ncols

    fig.text(0.008, y_of("title") - 0.02 / height, TITLE, fontsize=14.5, color=INK,
             ha="left", va="top")
    for i, line in enumerate(subtitle_lines):
        fig.text(0.008, y_of("subtitle") - (0.02 + 0.20 * i) / height, line,
                 fontsize=9.2, color=INK_SOFT, ha="left", va="top")

    fig.text(0.008, y_of("section_a") - 0.04 / height, "A · Runtime, split by phase",
             fontsize=10.5, color=INK, ha="left", va="top", fontweight="bold")
    for i, line in enumerate(note_a_lines):
        fig.text(0.008, y_of("section_a_note") - (0.02 + 0.19 * i) / height, line,
                 fontsize=8.6, color=INK_SOFT, ha="left", va="top")

    time_axes = []
    for i, (dataset, n_seqs) in enumerate(sizes):
        rect = [left_frac + i * (col_w + col_gap_in / width),
                block_bottom("panel_a") + axis_in / height,
                col_w, (panel_a_in - panel_title_in - axis_in) / height]
        ax = fig.add_axes(rect)
        draw_time_panel(ax, table, dataset, n_seqs, timeout_hint, show_labels=(i == 0))
        time_axes.append(ax)

    fig.legend(handles=[
        Patch(facecolor=RESIDUAL, label="other SATIVA steps (bookkeeping, rounding)"),
        Patch(facecolor=SHARED, label="reference tree (RAxML)"),
        Patch(facecolor=ORANGE, label="leave-one-out: RAxML -f O"),
        Patch(facecolor=AQUA, label="leave-one-out: EPA-ng"),
    ], loc="upper left", bbox_to_anchor=(left_frac, y_of("legend_a")), frameon=False,
        fontsize=8.5, ncol=legend_cols, labelcolor=INK, handlelength=1.1,
        handleheight=1.1, columnspacing=1.6, borderaxespad=0)

    fig.text(0.008, y_of("section_b") - 0.04 / height, "B · Per-sequence agreement with unmodified SATIVA",
             fontsize=10.5, color=INK, ha="left", va="top", fontweight="bold")
    # The control row is the point of the panel, and it needs saying where it is read.
    for i, line in enumerate(note_b_lines):
        fig.text(0.008, y_of("section_b_note") - (0.02 + 0.19 * i) / height, line,
                 fontsize=8.6, color=INK_SOFT, ha="left", va="top")
    ax_b = fig.add_axes([left_frac, block_bottom("panel_b") + axis_in / height,
                         x_right - left_frac,
                         (panel_b_in - 0.12 - axis_in) / height])
    draw_agreement_panel(ax_b, agreement_rows)

    fig.legend(handles=[
        Patch(facecolor=AQUA, label="same sequence, same rank"),
        Patch(facecolor=YELLOW, label="same sequence, other rank"),
        Patch(facecolor=SHARED, label="flagged by the original only"),
        Patch(facecolor=BLUE, label="flagged by EPA-ng only"),
    ], loc="upper left", bbox_to_anchor=(left_frac, y_of("legend_b")), frameon=False,
        fontsize=8.5, ncol=legend_cols, labelcolor=INK, handlelength=1.1,
        handleheight=1.1, columnspacing=1.6, borderaxespad=0)

    fig.text(0.008, y_of("section_c") - 0.04 / height,
             f"C · The same, keeping only calls above {strict_cutoff:g} confidence",
             fontsize=10.5, color=INK, ha="left", va="top", fontweight="bold")
    for i, line in enumerate(note_c_lines):
        fig.text(0.008, y_of("section_c_note") - (0.02 + 0.19 * i) / height, line,
                 fontsize=8.6, color=INK_SOFT, ha="left", va="top")
    ax_c = fig.add_axes([left_frac, block_bottom("panel_c") + axis_in / height,
                         x_right - left_frac,
                         (panel_c_in - 0.12 - axis_in) / height])
    draw_agreement_panel(ax_c, strict_agreement_rows)

    for i, line in enumerate(footer_lines):
        fig.text(0.008, y_of("footer") - (0.16 + 0.17 * i) / height, line,
                 fontsize=8, color=MUTED, ha="left", va="top")

    for suffix in ("png", "pdf"):
        path = outdir / f"fig_sativa_speedup.{suffix}"
        fig.savefig(path, dpi=200, facecolor="white")
        print(f"wrote {path}")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs-dir", default=str(BENCH_ROOT / "results" / "runs"))
    parser.add_argument("--outdir", default=str(BENCH_ROOT / "results"))
    parser.add_argument("--timeout-hint", type=float, default=7200)
    parser.add_argument("--min-size", type=int, default=MIN_SIZE_TIME_PANEL,
                        help="smallest alignment size shown in panel A")
    parser.add_argument("--strict-cutoff", type=float, default=0.9,
                        help="confidence filter applied to both sides in panel C")
    args = parser.parse_args()

    records = load_runs(args.runs_dir)
    if not records:
        raise SystemExit(f"no run.json under {args.runs_dir}")
    table = aggregate(records)

    outdir = Path(args.outdir)
    figures = outdir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    rows = write_summary(table, outdir / "summary.tsv")
    print(f"{len(records)} runs -> {len(rows)} cells -> {outdir / 'summary.tsv'}")

    figure_main(table,
                load_concordance(args.runs_dir, BASELINE),
                load_concordance(args.runs_dir, BASELINE, args.strict_cutoff),
                figures, args.timeout_hint, args.min_size, args.strict_cutoff)


if __name__ == "__main__":
    main()
