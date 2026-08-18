#!/usr/bin/env python3
# Part of 18.sativa_speedup. What the numbers mean: ../RESULTS.md
"""Build the size gradient for the SATIVA speed benchmark.

Draws nested random subsets from a pipeline MSA (headers
`>id|kingdom=..|...|species=..`) and writes, per size, the two files SATIVA
consumes: `<name>.fasta` (ids only) and `<name>.tax` (id<TAB>K;P;C;O;F;G;S).

Nested = the n=100 set is a subset of n=200, etc., so a runtime difference
between sizes comes from the number of sequences, not from a different draw.
Columns that are all-gap after subsetting are dropped: they carry no signal and
SATIVA/RAxML would reduce the alignment anyway (doing it here keeps the
alignment width honest and comparable across conditions).
"""

import argparse
import json
import random
import unicodedata
from pathlib import Path

TAX_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]
GAP = set("-.")

# The 2017 python-2 SATIVA aborts on any non-ASCII taxon name (UnicodeEncodeError in its
# bundled ete2 while parsing its own taxonomy Newick). In this ITS cluster the only
# offender is the U+00D7 of hybrid plant names -- "Ilex x altaclerensis" is the standard
# ASCII spelling of the same taxon -- so transliterating costs nothing and lets the
# unmodified python-2 build run the whole gradient as a control.
ASCII_MAP = {"\u00d7": "x", "\u00f7": "/"}


def read_msa(path):
    """Yield (header, sequence) pairs from a FASTA/aligned FASTA file."""
    header, parts = None, []
    with open(path) as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:], []
            else:
                parts.append(line)
    if header is not None:
        yield header, "".join(parts)


def to_ascii(value):
    """ASCII spelling of a taxon name; returns (text, changed)."""
    if all(ord(c) < 128 for c in value):
        return value, False
    text = "".join(ASCII_MAP.get(c, c) for c in value)
    text = "".join(c for c in unicodedata.normalize("NFKD", text)
                   if not unicodedata.combining(c))
    text = text.encode("ascii", "ignore").decode("ascii")
    text = " ".join(text.split())
    return text, True


def parse_header(header, ascii_taxonomy=False):
    """Split a pipeline MSA header into (seq_id, 7-rank taxonomy string)."""
    fields = header.split("|")
    seq_id = fields[0]
    info = {}
    for item in fields[1:]:
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
    # Same normalisation as workflow/scripts/Sativa/auto_sativa.py: underscores
    # were introduced to protect spaces through MAFFT, and an empty rank breaks
    # SATIVA's Newick taxonomy parsing.
    ranks, changed = [], False
    for rank in TAX_ORDER:
        value = (info.get(rank) or "").strip().replace("_", " ")
        if ascii_taxonomy:
            value, was_changed = to_ascii(value)
            changed = changed or was_changed
        ranks.append(value if value else "NA")
    return seq_id, ";".join(ranks), changed


def drop_all_gap_columns(seqs):
    """Remove the alignment columns that are gap-only in this subset."""
    width = len(seqs[0])
    keep = [i for i in range(width)
            if any(seq[i] not in GAP for seq in seqs)]
    if len(keep) == width:
        return seqs, width, width
    return ["".join(seq[i] for i in keep) for seq in seqs], width, len(keep)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", required=True, help="Source MSA (pipeline header format)")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--sizes", required=True,
                        help="Comma-separated subset sizes, e.g. 100,200,400,800,1600. "
                             "'full' adds the complete alignment.")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--ascii-taxonomy", action="store_true",
                        help="transliterate taxon names to ASCII (see ASCII_MAP): the "
                             "python-2 SATIVA cannot parse anything else")
    args = parser.parse_args()

    records = [(h, s) for h, s in read_msa(args.source)]
    total = len(records)
    print(f"source: {args.source}\n  {total} sequences x {len(records[0][1])} columns")

    rng = random.Random(args.seed)
    order = list(range(total))
    rng.shuffle(order)  # one shuffle -> every subset is a prefix -> nested subsets

    sizes = []
    for token in args.sizes.split(","):
        token = token.strip()
        if not token:
            continue
        sizes.append(total if token == "full" else int(token))
    sizes = sorted({s for s in sizes if s <= total})

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = []

    for size in sizes:
        name = f"n{size}"
        subset = [records[i] for i in order[:size]]
        seqs, width_in, width_out = drop_all_gap_columns([s for _, s in subset])

        dest = outdir / name
        dest.mkdir(exist_ok=True)
        fasta_path, tax_path = dest / f"{name}.fasta", dest / f"{name}.tax"
        n_transliterated = 0
        with open(fasta_path, "w") as fasta, open(tax_path, "w") as tax:
            for (header, _), seq in zip(subset, seqs):
                seq_id, taxonomy, changed = parse_header(header, args.ascii_taxonomy)
                n_transliterated += bool(changed)
                fasta.write(f">{seq_id}\n{seq}\n")
                tax.write(f"{seq_id}\t{taxonomy}\n")

        entry = {"name": name, "n_seqs": size, "aln_len": width_out,
                 "aln_len_before_gap_trim": width_in,
                 "ascii_taxonomy": bool(args.ascii_taxonomy),
                 "n_transliterated": n_transliterated,
                 "fasta": str(fasta_path), "tax": str(tax_path)}
        manifest.append(entry)
        extra = f", {n_transliterated} taxonomies transliterated" if n_transliterated else ""
        print(f"  {name}: {size} seqs x {width_out} cols (from {width_in}){extra}")

    meta = {"source": str(args.source), "source_n_seqs": total, "seed": args.seed,
            "nested": True, "ascii_taxonomy": bool(args.ascii_taxonomy),
            "datasets": manifest}
    (outdir / "datasets.json").write_text(json.dumps(meta, indent=2))
    print(f"wrote {outdir / 'datasets.json'}")


if __name__ == "__main__":
    main()
