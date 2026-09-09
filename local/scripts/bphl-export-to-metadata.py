#!/usr/bin/env python3
"""
Join the three BPHL exports into the metadata table the local workflow reads.

The sequencing results and the epidemiology arrive as separate tab-delimited
files sharing only a sample identifier:

    sequences.txt   sampleID, serotype, nextclade_clade
    metadata.txt    sampleID, Imported Status, Origin, Date of Collection,
                    Collection County
    mosquito.txt    sampleID, Species, Origin, Date of Collection

Output is defaults/metadata_template.tsv's schema, one row per sequenced sample
that has a collection date. Use summary-report-to-metadata.py instead when
converting a single Daytona_dengue run.
"""

import argparse
import csv
import re
import sys
from collections import OrderedDict

COLUMNS = [
    "sample_id",
    "serotype",
    "nextclade_clade",
    "collection_date",
    "location",
    "case_origin",
    "travel_country",
    "host",
    "strain",
    "notes",
]

# Origin holds a country for imported cases and "FL - <county>" for locally
# acquired ones, so the two have to be told apart before either is used.
FLORIDA = re.compile(r"^FL\b[\s-]*", re.IGNORECASE)
REGION_ONLY = {"South America", "Central America", "Africa", "Asia", "Caribbean"}
UNKNOWN_ORIGIN = {"", "Unknown", "unknown"}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--mosquito", required=True)
    parser.add_argument("--synonyms", required=True)
    parser.add_argument("--countries", help="color_orderings.tsv, to check travel_country")
    parser.add_argument("--output", required=True)
    parser.add_argument("--report", required=True)
    return parser.parse_args()


def read_table(path):
    # newline="" plus the csv module handles both the CRLF line endings and the
    # CSV-style quoting around the comma-bearing Imported Status values.
    with open(path, newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [
            {key: (value or "").strip() for key, value in row.items() if key}
            for row in reader
        ]


def read_synonyms(path):
    synonyms = {}
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if not line or line.startswith("#"):
                continue
            source, target = line.split("\t")
            synonyms[source.strip()] = target.strip()
    return synonyms


def read_countries(path):
    with open(path, encoding="utf-8") as handle:
        return {
            line.split("\t", 1)[1].strip()
            for line in handle
            if line.startswith("country\t")
        }


def join_keys(sample_id):
    """
    Candidate epi identifiers for a sequencing identifier, least stripped first.

    Re-sequenced samples carry decorations the epi export does not: a t_ prefix,
    a _NC_/_RJ_ run suffix, and a trailing K. Stripping is tried rather than
    applied so that an identifier genuinely ending in K is not mangled.
    """
    candidates = [sample_id]
    for candidate in list(candidates):
        if candidate.startswith("t_"):
            candidates.append(candidate[2:])
    for candidate in list(candidates):
        stripped = re.sub(r"_(NC|RJ)_\d+$", "", candidate)
        if stripped != candidate:
            candidates.append(stripped)
    for candidate in list(candidates):
        if candidate.endswith("K"):
            candidates.append(candidate[:-1])
    return list(OrderedDict.fromkeys(candidates))


def county_from_origin(origin):
    if not FLORIDA.match(origin):
        return ""
    county = FLORIDA.sub("", origin).strip()
    return "" if county.lower() in {"", "unknown county"} else county


def resolve_origin(origin, synonyms, countries, unresolved):
    """
    (travel_country, is_ambiguous) for an Origin value.

    Ambiguous means the case cannot be placed: several origins listed, a region
    rather than a country, or nothing recorded. Those become case_origin
    unknown rather than being resolved to whichever value happens to be first.
    """
    if origin in UNKNOWN_ORIGIN or ";" in origin or origin in REGION_ONLY:
        return "", True
    if FLORIDA.match(origin):
        return "", True

    country = synonyms.get(origin, origin)
    if countries is not None and country not in countries:
        unresolved.setdefault(country, 0)
        unresolved[country] += 1
    return country, False


def main():
    args = parse_args()

    synonyms = read_synonyms(args.synonyms)
    countries = read_countries(args.countries) if args.countries else None

    human = {row["sampleID"]: row for row in read_table(args.metadata)}
    mosquito = {row["sampleID"]: row for row in read_table(args.mosquito)}

    records = []
    undated = []
    unresolved = {}
    replicates = OrderedDict()

    for row in read_table(args.sequences):
        sample_id = row["sampleID"]
        epi = None
        source = None
        for key in join_keys(sample_id):
            if key in mosquito:
                epi, source, base = mosquito[key], "mosquito", key
                break
            if key in human:
                epi, source, base = human[key], "human", key
                break

        if epi is None:
            undated.append(sample_id)
            continue

        replicates.setdefault(base, []).append(sample_id)

        origin = epi.get("Origin", "")
        record = {column: "" for column in COLUMNS}
        record["sample_id"] = sample_id
        record["serotype"] = row["serotype"]
        record["nextclade_clade"] = row["nextclade_clade"]
        record["collection_date"] = epi["Date of Collection"]

        if source == "mosquito":
            record["host"] = epi.get("Species", "")
            record["location"] = county_from_origin(origin)
            record["case_origin"] = "local"
        else:
            record["location"] = epi.get("Collection County", "")
            status = epi.get("Imported Status", "")
            if status.startswith("Unknown"):
                record["case_origin"] = "unknown"
            elif status == "Acquired in Florida":
                record["case_origin"] = "local"
                exposure = county_from_origin(origin)
                if exposure and exposure != record["location"]:
                    record["notes"] = f"exposure county {exposure}"
            else:
                country, ambiguous = resolve_origin(
                    origin, synonyms, countries, unresolved
                )
                if ambiguous:
                    record["case_origin"] = "unknown"
                    record["notes"] = f"origin recorded as {origin!r}"
                else:
                    record["case_origin"] = "travel-associated"
                    record["travel_country"] = country

        records.append(record)

    with open(args.output, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=COLUMNS, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(records)

    groups = {base: ids for base, ids in replicates.items() if len(ids) > 1}

    with open(args.report, "w", encoding="utf-8") as handle:
        def emit(text):
            print(text, file=handle)
            print(text, file=sys.stderr)

        emit(f"wrote {len(records)} samples to {args.output}")
        emit(f"  travel-associated: {sum(1 for r in records if r['case_origin'] == 'travel-associated')}")
        emit(f"  local:             {sum(1 for r in records if r['case_origin'] == 'local')}")
        emit(f"  unknown:           {sum(1 for r in records if r['case_origin'] == 'unknown')}")

        if undated:
            emit(f"\ndropped {len(undated)} sequenced samples with no epidemiology row.")
            emit("They have no collection date, which is required. Supply one to include them.")
            for sample_id in sorted(undated):
                emit(f"  {sample_id}")

        if groups:
            emit(f"\n{len(groups)} samples were sequenced more than once. All rows are kept.")
            emit("Compare unambiguous base counts in results/validation_report.txt and")
            emit("delete the weaker row from the metadata before building.")
            for base, ids in groups.items():
                emit(f"  {base}: {', '.join(ids)}")

        if unresolved:
            emit(f"\n{len(unresolved)} travel countries are absent from the colour ordering file.")
            emit("They will fall outside every region group when colours are assigned.")
            for country, count in sorted(unresolved.items()):
                emit(f"  {country} ({count})")


if __name__ == "__main__":
    main()
