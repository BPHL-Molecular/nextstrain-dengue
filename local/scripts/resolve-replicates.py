#!/usr/bin/env python3
"""
Keep one sequence per specimen when a sample was sequenced more than once.

Replicates share an identifier once the t_ prefix, the run suffix and a trailing
K or k are removed. A VADR PASS copy is kept over any other; among copies with
the same flag, the one with the most unambiguous bases. The others are dropped
from both the metadata and the FASTA.
"""

import argparse
import csv
import re
import sys

from Bio import SeqIO

UNAMBIGUOUS = set("ACGTU")
RUN_SUFFIX = re.compile(r"(?:_(?:NC|RJ)_\d+|[-_](?:repeat\d*|NextSeq|test))$")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--id-column", default="accession")
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-sequences", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def undecorated(sample_id):
    # Keep in step with join_keys in prepare-daytona-inputs.py.
    stripped = sample_id[2:] if sample_id.startswith("t_") else sample_id
    return RUN_SUFFIX.sub("", stripped)


def main():
    args = parse_args()

    with open(args.metadata, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        columns = reader.fieldnames
        rows = list(reader)

    sequences = {record.id: str(record.seq) for record in SeqIO.parse(args.sequences, "fasta")}

    keys = {row[args.id_column]: undecorated(row[args.id_column]) for row in rows}
    flag = {row[args.id_column]: row.get("vadr_flag", "").upper() for row in rows}
    present = set(keys.values())
    # A trailing K only marks a replicate when the unsuffixed identifier exists too,
    # so an identifier that genuinely ends in K is left alone.
    for sample_id, key in keys.items():
        if key[-1:] in ("K", "k") and key[:-1] in present:
            keys[sample_id] = key[:-1]

    groups = {}
    for sample_id, key in keys.items():
        groups.setdefault(key, []).append(sample_id)

    def unambiguous(sample_id):
        return sum(1 for base in sequences.get(sample_id, "").upper() if base in UNAMBIGUOUS)

    dropped = set()
    report = []
    for key, members in groups.items():
        if len(members) < 2:
            continue
        ranked = sorted(members, key=lambda m: (flag[m] != "PASS", -unambiguous(m), len(m), m))
        dropped.update(ranked[1:])
        for rank, member in enumerate(ranked):
            report.append((key, member, unambiguous(member), "kept" if rank == 0 else "dropped"))

    with open(args.output_metadata, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(row for row in rows if row[args.id_column] not in dropped)

    with open(args.output_sequences, "w", encoding="utf-8") as out:
        for sample_id, sequence in sequences.items():
            if sample_id not in dropped:
                out.write(f">{sample_id}\n{sequence}\n")

    with open(args.report, "w", encoding="utf-8") as out:
        out.write("specimen\tsample_id\tunambiguous_bases\tstatus\n")
        for line in report:
            out.write("\t".join(map(str, line)) + "\n")

    print(f"{len(dropped)} replicate sequences dropped from {sum(len(m) > 1 for m in groups.values())} specimens", file=sys.stderr)


if __name__ == "__main__":
    main()
