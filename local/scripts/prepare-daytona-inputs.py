#!/usr/bin/env python3
"""
Prepare the local metadata and sequences from Daytona_dengue runs.

Each run folder holds summary_report.txt and assemblies_qc_pass/{pass,review}/
with one consensus FASTA per sample. Daytona_dengue writes every PASS assembly to
pass/ but only the REVIEW assemblies covering at least 80% of the genome to
review/, so the folders, not the flag alone, decide what is usable. A sample
enters the build when a run reports it with an accepted VADR flag, its FASTA is
in that flag's folder, and the sample metadata has a row for it.
"""

import argparse
import csv
import os
import re
import sys
from collections import Counter, defaultdict

from Bio import SeqIO

METADATA_COLUMNS = ["sample_id", "collection_date", "location", "case_origin", "travel_country"]
REPORT_COLUMNS = ["sample_id", "serotype", "nextclade_clade", "vadr_flag"]
OUTPUT_COLUMNS = [
    "sample_id",
    "serotype",
    "nextclade_clade",
    "collection_date",
    "location",
    "case_origin",
    "travel_country",
    "host",
    "vadr_flag",
]
RUN_SUFFIX = re.compile(r"(?:_(?:NC|RJ)_\d+|[-_](?:repeat\d*|NextSeq|test))$")
FASTA_NAME = re.compile(r"(?:[._]consensus)?\.(?:fasta|fa|fna)$", re.IGNORECASE)


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--runs", nargs="*", default=[], help="Daytona_dengue output folders, later ones win ties")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--vadr-flags", nargs="+", default=["PASS", "REVIEW"])
    parser.add_argument("--vector-pattern", default=r"(?i)mosquito")
    parser.add_argument("--vector-host", default="Aedes aegypti")
    parser.add_argument("--synonyms", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-sequences", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def read_tsv(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = [{(k or "").strip().lower(): (v or "").strip() for k, v in row.items() if k} for row in reader]
        header = [(name or "").strip().lower() for name in reader.fieldnames or []]
    return header, rows


def require_columns(path, header, columns):
    missing = [column for column in columns if column not in header]
    if missing:
        raise SystemExit(f"{path} is missing columns {missing}; found {header}")


def read_synonyms(path):
    synonyms = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line and not line.startswith("#"):
                source, target = line.split("\t")
                synonyms[source.strip()] = target.strip()
    return synonyms


def join_keys(sample_id):
    """Candidate metadata identifiers, least stripped first, for a sequencing run identifier."""
    candidates = [sample_id]
    for candidate in list(candidates):
        if candidate.startswith("t_"):
            candidates.append(candidate[2:])
    for candidate in list(candidates):
        stripped = RUN_SUFFIX.sub("", candidate)
        if stripped != candidate:
            candidates.append(stripped)
    for candidate in list(candidates):
        if candidate[-1:] in ("K", "k"):
            candidates.append(candidate[:-1])
    return list(dict.fromkeys(candidates))


def index_fasta(folder):
    """Every FASTA in folder, keyed by its file name stem and by its first header."""
    index = defaultdict(set)
    if not os.path.isdir(folder):
        return index
    for name in sorted(os.listdir(folder)):
        if not FASTA_NAME.search(name):
            continue
        path = os.path.join(folder, name)
        index[FASTA_NAME.sub("", name)].add(path)
        first = next(SeqIO.parse(path, "fasta"), None)
        if first is not None:
            index[first.id].add(path)
    return index


def main():
    args = parse_args()
    if not args.runs:
        raise SystemExit("no Daytona runs given; list their output folders under daytona_runs")

    accepted = {flag.upper() for flag in args.vadr_flags}
    vector = re.compile(args.vector_pattern)
    synonyms = read_synonyms(args.synonyms)

    header, metadata_rows = read_tsv(args.metadata)
    require_columns(args.metadata, header, METADATA_COLUMNS)
    metadata = {}
    for row in metadata_rows:
        if row["sample_id"] in metadata:
            raise SystemExit(f"{args.metadata} lists {row['sample_id']!r} more than once")
        if row["sample_id"]:
            metadata[row["sample_id"]] = row

    flag_counts = Counter()
    chosen = {}
    flags_by_metadata_row = defaultdict(set)
    no_serotype, missing_fasta, not_in_folder, several_files, duplicates = [], [], [], [], []

    for order, run in enumerate(args.runs):
        report_path = os.path.join(run, "summary_report.txt")
        header, rows = read_tsv(report_path)
        require_columns(report_path, header, REPORT_COLUMNS)
        folders = {
            "PASS": index_fasta(os.path.join(run, "assemblies_qc_pass", "pass")),
            "REVIEW": index_fasta(os.path.join(run, "assemblies_qc_pass", "review")),
        }

        for row in rows:
            sample_id, flag = row["sample_id"], row["vadr_flag"].upper()
            flag_counts[flag] += 1
            for key in join_keys(sample_id):
                if key in metadata:
                    flags_by_metadata_row[key].add(flag)
                    break
            if flag not in accepted:
                continue
            if not re.search(r"[1-4]", row["serotype"]):
                no_serotype.append((run, sample_id, row["serotype"] or "NA"))
                continue

            paths = folders.get(flag, {}).get(sample_id, set())
            if not paths:
                if flag == "PASS":
                    missing_fasta.append((run, sample_id, flag))
                else:
                    not_in_folder.append((run, sample_id, flag))
                continue
            if len(paths) > 1:
                several_files.append((run, sample_id, sorted(paths)))
                continue

            candidate = {
                "run": run,
                "flag": flag,
                "serotype": row["serotype"],
                "nextclade_clade": row["nextclade_clade"],
                "path": next(iter(paths)),
            }
            previous = chosen.get(sample_id)
            if previous:
                winner = candidate if (flag == "PASS") >= (previous["flag"] == "PASS") else previous
                duplicates.append((sample_id, previous, candidate, winner))
                chosen[sample_id] = winner
            else:
                chosen[sample_id] = candidate

    kept, no_metadata, multi_record, used = [], [], [], set()
    with open(args.output_sequences, "w", encoding="utf-8") as fasta:
        for sample_id, run_row in chosen.items():
            key = next((k for k in join_keys(sample_id) if k in metadata), None)
            if key is None:
                no_metadata.append((sample_id, run_row["flag"]))
                continue
            records = list(SeqIO.parse(run_row["path"], "fasta"))
            if len(records) != 1:
                multi_record.append((sample_id, run_row["path"], len(records)))
                continue
            used.add(key)
            epi = metadata[key]
            travel_country = epi["travel_country"]
            kept.append({
                "sample_id": sample_id,
                "serotype": run_row["serotype"],
                "nextclade_clade": run_row["nextclade_clade"],
                "collection_date": epi["collection_date"],
                "location": epi["location"],
                "case_origin": epi["case_origin"],
                "travel_country": synonyms.get(travel_country, travel_country),
                "host": args.vector_host if vector.search(sample_id) else "",
                "vadr_flag": run_row["flag"],
            })
            fasta.write(f">{sample_id}\n{records[0].seq}\n")

    with open(args.output_metadata, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kept)

    ignored = [
        (sample_id, "not in any run" if sample_id not in flags_by_metadata_row
         else "VADR " + "/".join(sorted(flags_by_metadata_row[sample_id])))
        for sample_id in metadata if sample_id not in used
    ]

    with open(args.report, "w", encoding="utf-8") as handle:
        def emit(text=""):
            print(text, file=handle)
            print(text, file=sys.stderr)

        def section(title, items, show):
            if items:
                emit(f"\n{title}: {len(items)}")
                for item in items:
                    emit("  " + show(item))

        emit(f"runs read: {len(args.runs)}")
        for run in args.runs:
            emit(f"  {run}")
        emit("samples in the summary reports by VADR flag: "
             + ", ".join(f"{flag} {count}" for flag, count in sorted(flag_counts.items())))
        emit(f"accepted flags: {', '.join(sorted(accepted))}")
        emit(f"samples written: {len(kept)}"
             f" ({', '.join(f'{f} {n}' for f, n in sorted(Counter(r['vadr_flag'] for r in kept).items()))})")
        emit(f"metadata rows: {len(metadata)}, used: {len(used)}")

        section("excluded, no metadata row (no collection date)", no_metadata, lambda i: f"{i[0]}\t{i[1]}")
        section("excluded, PASS sample with no FASTA in pass/", missing_fasta, lambda i: f"{i[1]}\t{i[2]}\t{i[0]}")
        section("excluded, not in its flag's folder (REVIEW below the 80% coverage cutoff)", not_in_folder,
                lambda i: f"{i[1]}\t{i[2]}\t{i[0]}")
        section("excluded, more than one FASTA file matches", several_files, lambda i: f"{i[1]}\t{i[0]}\t{', '.join(i[2])}")
        section("excluded, FASTA file does not hold exactly one sequence", multi_record, lambda i: f"{i[0]}\t{i[1]}\t{i[2]} records")
        section("excluded, no assignable serotype", no_serotype, lambda i: f"{i[1]}\t{i[2]}\t{i[0]}")
        section("sequenced in more than one run", duplicates,
                lambda i: f"{i[0]}\tkept {i[3]['flag']} from {i[3]['run']}; "
                          f"runs {i[1]['flag']} {i[1]['run']} and {i[2]['flag']} {i[2]['run']}")
        section("metadata rows not used", ignored, lambda i: f"{i[0]}\t{i[1]}")

    if not kept:
        raise SystemExit("no samples passed every check; read the report above")


if __name__ == "__main__":
    main()
