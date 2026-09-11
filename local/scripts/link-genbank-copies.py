#!/usr/bin/env python3
"""
Link local samples to the copies of them already deposited in GenBank.

A BPHL genome submitted to NCBI returns through ingest under a GenBank accession,
so without this step the build carries it twice. Each alignment holds one
serotype's local sequences together with the GenBank candidates. A local
sequence and a GenBank record are the same genome when they agree at every
position where both have a called base and those positions cover most of the
shorter sequence.

A linked local row takes the GenBank accession and loses its FASTA record, so
augur merge joins it onto the GenBank row: the tip keeps the GenBank sequence and
accession while the local values, travel country and case origin among them,
overwrite the public ones.
"""

import argparse
import csv
import sys

import numpy as np
from Bio import SeqIO

CODES = np.zeros(256, dtype=np.uint8)
for code, base in enumerate("ACGT", start=1):
    CODES[ord(base)] = CODES[ord(base.lower())] = code


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--public-metadata", required=True)
    parser.add_argument("--alignments", nargs="+", required=True)
    parser.add_argument("--min-overlap", type=float, default=0.9)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-sequences", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def read_tsv(path):
    with open(path, newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def encode(sequence):
    return CODES[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]


def identical_pairs(alignment, local_ids, min_overlap):
    local, public = {}, {}
    for record in SeqIO.parse(alignment, "fasta"):
        (local if record.id in local_ids else public)[record.id] = encode(str(record.seq))
    if not local or not public:
        return []

    public_ids = list(public)
    matrix = np.vstack([public[i] for i in public_ids])
    called = matrix > 0
    public_called = called.sum(axis=1)

    pairs = []
    for local_id, seq in local.items():
        both = called & (seq > 0)
        shared = both.sum(axis=1)
        mismatches = (both & (matrix != seq)).sum(axis=1)
        needed = min_overlap * np.minimum(public_called, (seq > 0).sum())
        for index in np.flatnonzero((mismatches == 0) & (shared > 0) & (shared >= needed)):
            pairs.append((local_id, public_ids[index], int(shared[index])))
    return pairs


def narrow(options, date, date_of):
    same_day = [option for option in options if date_of[option] == date]
    return same_day if len(same_day) == 1 else options


def main():
    args = parse_args()

    columns, rows = read_tsv(args.metadata)
    local_by_id = {row["accession"]: row for row in rows}
    _, public_rows = read_tsv(args.public_metadata)
    public_by_id = {row["accession"]: row for row in public_rows}

    pairs = []
    for alignment in args.alignments:
        pairs.extend(identical_pairs(alignment, set(local_by_id), args.min_overlap))

    local_date = {i: row.get("date", "") for i, row in local_by_id.items()}
    public_date = {i: row.get("date", "") for i, row in public_by_id.items()}
    by_local, by_public = {}, {}
    for local_id, accession, _ in pairs:
        by_local.setdefault(local_id, []).append(accession)
        by_public.setdefault(accession, []).append(local_id)
    by_local = {i: narrow(o, local_date[i], public_date) for i, o in by_local.items()}
    by_public = {a: narrow(o, public_date.get(a, ""), local_date) for a, o in by_public.items()}

    links = {
        local_id: accession
        for local_id, accession, _ in pairs
        if by_local[local_id] == [accession] and by_public[accession] == [local_id]
    }
    linked_accessions = set(links.values())

    report = []
    for local_id, accession, shared in pairs:
        if links.get(local_id) == accession:
            status = "linked"
        elif local_id in links and accession in linked_accessions:
            continue
        else:
            status = "ambiguous"
        public = public_by_id.get(accession, {})
        report.append((status, local_id, accession, public.get("strain", ""),
                       local_by_id[local_id].get("serotype_genbank", ""), shared,
                       local_date[local_id], public_date.get(accession, "")))

    for accession, public in public_by_id.items():
        if "FL-BPHL" in public.get("strain", "") and accession not in linked_accessions:
            report.append(("bphl_named_unmatched", "", accession, public["strain"],
                           public.get("serotype_genbank", ""), "", "", public.get("date", "")))

    for local_id, accession in links.items():
        row = local_by_id[local_id]
        row["accession"] = accession
        row["accession_version"] = public_by_id.get(accession, {}).get("accession_version") or accession

    with open(args.output_metadata, "w", newline="", encoding="utf-8") as out:
        writer = csv.DictWriter(out, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    with open(args.output_sequences, "w", encoding="utf-8") as out:
        for record in SeqIO.parse(args.sequences, "fasta"):
            if record.id not in links:
                out.write(f">{record.id}\n{record.seq}\n")

    with open(args.report, "w", encoding="utf-8") as out:
        out.write("status\tlocal_id\taccession\tgenbank_strain\tserotype\tshared_bases\tlocal_date\tgenbank_date\n")
        for line in sorted(report, key=lambda r: (r[0], r[4], r[1], r[2])):
            out.write("\t".join(map(str, line)) + "\n")

    ambiguous = sum(1 for line in report if line[0] == "ambiguous")
    print(f"{len(links)} local samples linked to their GenBank copies", file=sys.stderr)
    if ambiguous:
        print(
            f"WARNING: {ambiguous} ambiguous pairs left unlinked, so both copies stay in "
            f"the build. Decide each row in {args.report} before running phylogenetic/.",
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()
