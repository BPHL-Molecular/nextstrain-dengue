#!/usr/bin/env python3
"""
Validate the hand-filled local metadata TSV against the consensus FASTA and emit
a metadata table with the same columns ingest produces.
"""

import argparse
import csv
import json
import os
import re
import sys
from datetime import date, datetime

from Bio import SeqIO

UNAMBIGUOUS = set("ACGTU")
VALID_NT = set("ACGTURYSWKMBDHVN-.")
CASE_ORIGINS = {"local", "travel-associated", "undetermined", ""}


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--sequences", required=True)
    parser.add_argument("--constants", required=True, help="JSON object of constant column values")
    parser.add_argument("--rules", required=True, help="JSON object of validation rules")
    parser.add_argument("--strain-template", required=True)
    parser.add_argument("--output-metadata", required=True)
    parser.add_argument("--output-report", required=True)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


def read_metadata(path):
    with open(path, newline="", encoding="utf-8-sig") as handle:
        lines = [line for line in handle if not line.lstrip().startswith("#")]
    reader = csv.DictReader(lines, delimiter="\t")
    if reader.fieldnames is None:
        raise SystemExit(f"{path} is empty")
    columns = [name.strip() for name in reader.fieldnames]
    rows = []
    for row in reader:
        rows.append({key.strip(): (value or "").strip() for key, value in row.items() if key is not None})
    return columns, rows


def read_sequences(path):
    records = {}
    duplicates = []
    for record in SeqIO.parse(path, "fasta"):
        seq_id = record.id.split()[0].strip()
        if seq_id in records:
            duplicates.append(seq_id)
        records[seq_id] = str(record.seq).upper()
    return records, duplicates


def read_public_accessions(path):
    """Accessions already in the public data, if ingest has been run."""
    if not path or not os.path.exists(path):
        return set()
    with open(path, newline="", encoding="utf-8") as handle:
        return {row["accession"] for row in csv.DictReader(handle, delimiter="\t")}


def read_host_map(path):
    """host -> (host_genus, host_type), from the same TSV ingest curates with."""
    if not path or not os.path.exists(path):
        return {}
    with open(path, newline="", encoding="utf-8") as handle:
        return {
            row["host"]: (row["host_genus"], row["host_type"])
            for row in csv.DictReader(handle, delimiter="\t")
        }


def read_region_map(path):
    """
    country -> region, read from the Auspice colour ordering file, where each
    "# Region" comment heads the countries belonging to it. Using that file
    rather than a second table keeps exposure regions on the same scale the
    colourings already use.
    """
    if not path or not os.path.exists(path):
        return {}
    mapping = {}
    region = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                region = line.lstrip("# ").strip()
            elif line.startswith("country\t") and region:
                mapping[line.split("\t", 1)[1].strip()] = region
    return mapping


def normalise_serotype(value):
    match = re.search(r"([1-4])", value or "")
    if not match:
        return ""
    return f"denv{match.group(1)}"


def parse_date(value):
    """Accept a full ISO date or a partial date with XX for unknown parts."""
    if re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        return datetime.strptime(value, "%Y-%m-%d").date(), int(value[:4])
    if re.fullmatch(r"\d{4}-\d{2}-XX", value):
        return None, int(value[:4])
    if re.fullmatch(r"\d{4}-XX-XX", value):
        return None, int(value[:4])
    raise ValueError(value)


def sequence_stats(sequence):
    ungapped = sequence.replace("-", "").replace(".", "")
    unambiguous = sum(1 for base in ungapped if base in UNAMBIGUOUS)
    ambiguous_fraction = 0.0 if not ungapped else 1 - (unambiguous / len(ungapped))
    return len(ungapped), unambiguous, ambiguous_fraction


def main():
    args = parse_args()
    constants = json.loads(args.constants)
    rules = json.loads(args.rules)
    strict = args.strict or rules.get("strict", False)

    errors = []
    warnings = []
    excluded = []

    columns, rows = read_metadata(args.metadata)
    sequences, duplicate_headers = read_sequences(args.sequences)

    known = set(rules["required_columns"]) | set(rules["optional_columns"])
    for column in columns:
        if column not in known:
            errors.append(f"unknown column {column!r}; expected one of {sorted(known)}")
    for column in rules["required_columns"]:
        if column not in columns:
            errors.append(f"required column {column!r} is missing")

    for seq_id in sorted(set(duplicate_headers)):
        errors.append(f"duplicate FASTA record {seq_id!r}")

    if errors:
        report_and_exit(args, errors, warnings, [])

    id_pattern = re.compile(rules["id_regex"])
    host_map = read_host_map(rules.get("host_map"))
    region_map = read_region_map(rules.get("region_map"))
    allowed = set(rules["allowed_serotypes"])
    today = date.today()

    seen = set()
    records = []
    for index, row in enumerate(rows, start=2):
        where = f"row {index}"
        sample_id = row.get("sample_id", "")

        if not sample_id:
            errors.append(f"{where}: empty sample_id")
            continue
        if sample_id in seen:
            errors.append(f"{where}: duplicate sample_id {sample_id!r}")
            continue
        seen.add(sample_id)
        if not id_pattern.fullmatch(sample_id):
            errors.append(f"{where}: sample_id {sample_id!r} does not match {rules['id_regex']}")
            continue

        serotype = normalise_serotype(row.get("serotype", ""))
        if serotype not in allowed:
            errors.append(f"{sample_id}: serotype {row.get('serotype', '')!r} is not one of {sorted(allowed)}")
            continue

        raw_date = row.get("collection_date", "")
        try:
            parsed, year = parse_date(raw_date)
        except ValueError:
            errors.append(
                f"{sample_id}: collection_date {raw_date!r} is not YYYY-MM-DD, YYYY-MM-XX, or YYYY-XX-XX"
            )
            continue
        if parsed and parsed > today:
            errors.append(f"{sample_id}: collection_date {raw_date} is in the future")
            continue
        if year < 1950:
            errors.append(f"{sample_id}: collection_date {raw_date} predates 1950")
            continue

        record = dict(constants)
        record.update({key: value for key, value in row.items() if value})

        # augur filter drops rows where any of these is a literal "?".
        for field in ("country", "region"):
            if record.get(field, "") in ("", "?"):
                errors.append(f"{sample_id}: {field} resolved to {record.get(field, '')!r}")

        # Mirrors ingest's transform-new-fields --pass-through true: an unmapped
        # host becomes its own genus and type rather than being dropped.
        host = record.get("host", "")
        host_genus, host_type = host_map.get(host, (host, host))
        if host not in host_map:
            warnings.append(
                f"{sample_id}: host {host!r} is not in the host map, so host_genus "
                "and host_type fall back to the host name"
            )

        # Daytona already called this against the v-gen-lab datasets, in the same
        # format ingest writes, so it is carried through rather than recomputed.
        clade = row.get("nextclade_clade", "")
        if not clade or not re.match(r"^[0-9][A-Z]", clade):
            warnings.append(
                f"{sample_id}: nextclade_clade {clade or 'empty'!r} is not a dengue "
                "lineage, so the lineage colorings will be blank for this sample"
            )
            clade = ""

        # augur traits reconstructs ancestral geography from the exposure
        # columns, so an imported case has to point at where infection happened
        # rather than where it was reported, or Florida reads as a source.
        travel_country = record.get("travel_country", "")
        country_exposure = travel_country or record.get("country", "")
        region_exposure = record.get("region", "")
        if travel_country:
            region_exposure = region_map.get(travel_country, region_exposure)
            if travel_country not in region_map:
                warnings.append(
                    f"{sample_id}: travel_country {travel_country!r} is not in the colour "
                    f"ordering file, so region_exposure stays {region_exposure!r}"
                )

        record.update(
            {
                "accession": sample_id,
                "accession_version": sample_id,
                "date": raw_date,
                "serotype_genbank": serotype,
                "host": host,
                "host_genus": host_genus,
                "host_type": host_type,
                "genotype_nextclade": clade,
                "country_exposure": country_exposure,
                "region_exposure": region_exposure,
                "date_released": today.isoformat(),
                "date_updated": today.isoformat(),
            }
        )
        for consumed in ("sample_id", "collection_date", "serotype", "nextclade_clade"):
            record.pop(consumed, None)

        if not record.get("strain"):
            record["strain"] = args.strain_template.format(
                serotype_number=serotype[-1],
                country=record["country"].replace(" ", "_").upper(),
                accession=sample_id,
                year=year,
            )

        case_origin = record.get("case_origin", "")
        if case_origin and case_origin not in CASE_ORIGINS:
            warnings.append(
                f"{sample_id}: case_origin {case_origin!r} is not one of {sorted(CASE_ORIGINS)}"
            )

        if travel_country:
            if not case_origin:
                # Recording a travel country is what travel-associated means.
                case_origin = "travel-associated"
            elif case_origin != "travel-associated":
                warnings.append(
                    f"{sample_id}: travel_country is set but case_origin is {case_origin!r}"
                )
            if host_type and host_type != "Human":
                warnings.append(
                    f"{sample_id}: travel_country is set on a {host_type.lower()} sample. "
                    "Vectors are collected where they are found, so the exposure country "
                    "should be the collection country and a suspected origin belongs in notes"
                )
        elif not case_origin:
            # Deliberately not defaulted to "local". An unrecorded travel history
            # and an uninvestigated one look identical here, and asserting local
            # transmission that was never established is the worse error.
            warnings.append(
                f"{sample_id}: case_origin is blank, so it will render as an empty "
                "category in Auspice. Set it to 'local' or 'undetermined'"
            )

        record["case_origin"] = case_origin

        records.append(record)

    metadata_ids = {record["accession"] for record in records}
    for sample_id in sorted(metadata_ids - set(sequences)):
        errors.append(f"{sample_id}: present in metadata but has no FASTA record")

    # The FASTA is allowed to be a superset. Concatenating every consensus and
    # then listing only the QC-passing ones in the metadata is the normal way to
    # select what enters the trees, so these are recorded rather than warned
    # about: under strict they must not fail the run.
    excluded.extend(sorted(set(sequences) - metadata_ids))
    if excluded:
        print(f"Excluding {len(excluded)} FASTA records absent from the metadata.", file=sys.stderr)

    for sample_id in sorted(metadata_ids & read_public_accessions(rules.get("public_metadata"))):
        errors.append(
            f"{sample_id}: collides with a public GenBank accession; "
            "augur merge would silently keep only one of the two sequences"
        )

    stats = []
    for record in records:
        sequence = sequences.get(record["accession"])
        if sequence is None:
            continue
        invalid = set(sequence) - VALID_NT
        if invalid:
            errors.append(f"{record['accession']}: sequence contains {sorted(invalid)}")
        length, unambiguous, ambiguous_fraction = sequence_stats(sequence)
        record["length"] = str(length)
        stats.append((record["accession"], record["serotype_genbank"], length, unambiguous, ambiguous_fraction))

        if unambiguous < rules["min_ungapped_length"]:
            warnings.append(
                f"{record['accession']}: {unambiguous} unambiguous bases is below "
                f"{rules['min_ungapped_length']}; it will only enter the genome tree because "
                "include.txt forces it past --min-length"
            )
        if unambiguous < 1000:
            warnings.append(
                f"{record['accession']}: {unambiguous} unambiguous bases is below 1000; "
                "E gene extraction will drop it"
            )
        if ambiguous_fraction > rules["max_ambiguous_fraction"]:
            warnings.append(
                f"{record['accession']}: {ambiguous_fraction:.1%} ambiguous exceeds "
                f"{rules['max_ambiguous_fraction']:.0%}"
            )

    if errors or (strict and warnings):
        report_and_exit(args, errors, warnings, stats, excluded)

    write_report(args.output_report, errors, warnings, stats, excluded)

    fieldnames = sorted({key for record in records for key in record})
    with open(args.output_metadata, "w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for record in records:
            writer.writerow(record)

    print(f"Validated {len(records)} local samples.", file=sys.stderr)


def write_report(path, errors, warnings, stats, excluded=()):
    lines = ["# Local metadata validation report", ""]

    by_serotype = {}
    for _, serotype, _, _, _ in stats:
        by_serotype[serotype] = by_serotype.get(serotype, 0) + 1
    lines.append(f"Samples: {len(stats)}")
    for serotype in sorted(by_serotype):
        lines.append(f"  {serotype}: {by_serotype[serotype]}")
    lines.append("")

    lines.append("sample\tserotype\tlength\tunambiguous\tambiguous_fraction")
    for accession, serotype, length, unambiguous, fraction in sorted(stats):
        lines.append(f"{accession}\t{serotype}\t{length}\t{unambiguous}\t{fraction:.4f}")
    lines.append("")

    lines.append(f"Errors: {len(errors)}")
    lines.extend(f"  ERROR {message}" for message in errors)
    lines.append("")
    lines.append(f"Warnings: {len(warnings)}")
    lines.extend(f"  WARNING {message}" for message in warnings)
    lines.append("")
    lines.append(f"Excluded, present in the FASTA but not the metadata: {len(excluded)}")
    lines.extend(f"  {seq_id}" for seq_id in excluded)
    lines.append("")

    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def report_and_exit(args, errors, warnings, stats, excluded=()):
    write_report(args.output_report, errors, warnings, stats, excluded)
    for message in errors:
        print(f"ERROR {message}", file=sys.stderr)
    for message in warnings:
        print(f"WARNING {message}", file=sys.stderr)
    print(f"\nValidation failed. Full report: {args.output_report}", file=sys.stderr)
    sys.exit(1)


if __name__ == "__main__":
    main()
