#!/usr/bin/env python3
"""
Assign colours to high-cardinality metadata columns.

assign-colors.py gives a trait one ramp sized to its total number of values.
That works for a handful of regions or serotypes, but a build carrying 128
countries gets a 128-point interpolation in which neighbours differ by a couple
of units in one channel, and the ordering file puts the Americas at the end so
they land in a narrow orange-to-red band.

This splits countries by region first and gives each region its own full ramp, so
a region's members are spread across the whole spectrum rather than a sliver of
it. Colours then repeat between regions, which is the deliberate trade: region is
a separate colouring, and telling Cuba from Colombia matters more here than
telling Cuba from Cambodia.

The grouped columns share one assignment built from the union of their values, so
a country has the same colour under country and country_exposure.

Columns with no grouping to exploit, such as location, get a single ramp with the
same colour cycling assign-colors.py uses when values outnumber the palette.
"""

import argparse
import csv
import sys
from collections import OrderedDict

UNGROUPED = "Other"


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--ordering", required=True, help="color_orderings.tsv")
    parser.add_argument("--color-schemes", required=True)
    parser.add_argument(
        "--grouped-columns",
        nargs="+",
        default=["country", "country_exposure"],
        help="columns whose values are grouped by region before being coloured",
    )
    parser.add_argument(
        "--flat-columns",
        nargs="+",
        default=["location"],
        help="columns coloured from a single ramp, with cycling if needed",
    )
    parser.add_argument("--output", default="-")
    return parser.parse_args()


def read_region_of_country(path):
    """
    country -> region, from the "# Region" comment headers in the ordering file.

    Also returns the ordering position of each country so that values keep the
    file's geographic sequence within a region rather than falling into
    alphabetical order.
    """
    with open(path, encoding="utf-8") as handle:
        lines = handle.read().splitlines()

    # A comment counts as a region header only if it names one of the regions
    # the file itself declares. Section separators, prose comments and the
    # "# Hill 2024 dataset" heading are then ignored rather than being mistaken
    # for a region and mislabelling whatever countries follow them.
    regions = {
        line.split("\t", 1)[1].strip() for line in lines if line.startswith("region\t")
    }

    region_of = {}
    position = {}
    region = UNGROUPED
    for index, line in enumerate(lines):
        if line.startswith("#"):
            label = line.lstrip("#").strip()
            if label in regions:
                region = label
        elif line.startswith("country\t"):
            country = line.split("\t", 1)[1].strip()
            region_of.setdefault(country, region)
            position.setdefault(country, index)
    return region_of, position


def read_schemes(path):
    """Line N of the schemes file is an N-colour palette."""
    schemes = {}
    with open(path, encoding="utf-8") as handle:
        for number, line in enumerate(handle, start=1):
            line = line.strip()
            # A trailing blank line would otherwise register as a palette of one
            # empty string and become the largest scheme.
            if line:
                schemes[number] = line.split("\t")
    return schemes


def palette(schemes, count):
    """A palette of `count` colours, cycling the largest scheme if need be."""
    if count < 1:
        return []
    largest = max(schemes)
    if count <= largest:
        return schemes[count]

    colors = []
    remaining = count
    while remaining > largest:
        colors.extend(schemes[largest])
        remaining -= largest
    colors.extend(schemes[remaining])
    return colors


def observed_values(rows, column):
    """Distinct non-empty values, in first-seen order."""
    values = OrderedDict()
    for row in rows:
        value = (row.get(column) or "").strip()
        if value and value != "?":
            values[value] = None
    return list(values)


def main():
    args = parse_args()

    region_of, position = read_region_of_country(args.ordering)
    schemes = read_schemes(args.color_schemes)

    # csv handles the field quoting augur merge writes, so "Sri Lanka" arrives
    # as Sri Lanka and matches what augur export looks up.
    with open(args.metadata, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    out = sys.stdout if args.output == "-" else open(args.output, "w", encoding="utf-8")
    try:
        groups = OrderedDict()
        for column in args.grouped_columns:
            for value in observed_values(rows, column):
                groups.setdefault(region_of.get(value, UNGROUPED), OrderedDict())[value] = None

        color_of = {}
        for region, values in groups.items():
            ordered = sorted(values, key=lambda v: (position.get(v, len(position)), v))
            color_of.update(zip(ordered, palette(schemes, len(ordered))))
            print(f"{len(ordered)} values in {region}", file=sys.stderr)

        for column in args.grouped_columns:
            values = observed_values(rows, column)
            if not values:
                print(f"{column}: no values found, skipping", file=sys.stderr)
            for value in values:
                out.write(f"{column}\t{value}\t{color_of[value]}\n")

        for column in args.flat_columns:
            values = sorted(observed_values(rows, column))
            if not values:
                print(f"{column}: no values found, skipping", file=sys.stderr)
                continue
            for value, color in zip(values, palette(schemes, len(values))):
                out.write(f"{column}\t{value}\t{color}\n")
            print(f"{column}: {len(values)} values, ungrouped", file=sys.stderr)
    finally:
        if out is not sys.stdout:
            out.close()


if __name__ == "__main__":
    main()
