# local

Prepares privately sequenced dengue consensus genomes so the phylogenetic
workflow can place them on a tree alongside public GenBank data.

This workflow validates your metadata table, derives the dengue lineage levels
from the Nextclade call `Daytona_dengue` already made, and writes one metadata
file, one sequence file, and one include list per serotype.

## What you provide

Two files, both in `local/input/`, which is gitignored so consensus sequences
cannot be committed by accident.

### `local/input/sequences.fasta`

Every consensus genome for this run, concatenated into one FASTA. From a
`Daytona_dengue` run that is:

```sh
cat /path/to/daytona_output/assemblies_qc_pass/*.fasta > local/input/sequences.fasta
```

Each record header must be the sample identifier, either alone or as the first
whitespace-delimited token. The workflow rewrites headers to the bare
identifier, matching the bare accessions `ingest` writes.

### `local/input/metadata.tsv`

Generate it from the `Daytona_dengue` run rather than typing it out:

```sh
python3 local/scripts/summary-report-to-metadata.py \
    --summary-report /path/to/daytona_output/summary_report.txt \
    --output local/input/metadata.tsv
```

That fills in `sample_id`, `serotype`, and `nextclade_clade`, all of which the
pipeline already determined, and leaves `collection_date` and the epidemiological
columns for you. It keeps only samples whose `vadr_flag` is `PASS` and prints every sample it dropped, with
the reason. To include the `REVIEW` tier as well:

```sh
    --vadr-flags PASS,REVIEW
```

Samples whose serotype is `unclassified` are always dropped, because there is no
v-gen-lab dataset to place them against.

Mosquito pools are detected from the sample identifier and get `host` set to
`Aedes aegypti` instead of the human default, which the workflow then resolves to
`Aedes` and `Mosquito` through ingest's own host map. The converter names every
sample it treated this way. Adjust with `--vector-pattern` and `--vector-host` if
your identifiers or species differ; the host value has to appear in
`ingest/defaults/host_hostgenus_hosttype_map.tsv` or you get a warning and the
genus and type fall back to the host name.

The FASTA may contain more samples than the metadata. Anything present in the
sequences but absent from the metadata is excluded and listed in the validation
report, so concatenating all of `assemblies_qc_pass/` and then letting this
script decide what enters the trees is the intended workflow.

#### From a multi-year epidemiology export

When the sequencing results and the epidemiology come from separate exports
rather than a single `Daytona_dengue` run, use the other converter:

```sh
python3 local/scripts/bphl-export-to-metadata.py \
    --sequences data/sequences.txt \
    --metadata data/metadata.txt \
    --mosquito data/mosquito.txt \
    --synonyms local/defaults/country_synonyms.tsv \
    --countries phylogenetic/defaults/color_orderings.tsv \
    --output local/input/metadata.tsv \
    --report local/input/conversion_report.txt
```

It expects `sampleID`, `serotype` and `nextclade_clade` in the sequencing file;
`sampleID`, `Imported Status`, `Origin`, `Date of Collection` and
`Collection County` in the case file; and `sampleID`, `Species`, `Origin` and
`Date of Collection` in the mosquito file. Membership of the mosquito file is
what sets `host`, so a vector sample keeps its normal laboratory identifier
rather than needing a recognizable prefix.

`Imported Status` becomes `case_origin`, and `Origin` becomes `travel_country`
after passing through `defaults/country_synonyms.tsv`, which reconciles the
export's spellings with `phylogenetic/defaults/color_orderings.tsv`. An `Origin`
that names several places, a region rather than a country, or nothing at all
yields `case_origin: undetermined` and a blank `travel_country`, with the original
string kept in `notes`. For a locally acquired case `Origin` is the county of
exposure; it is recorded in `notes` when it differs from the collection county.

Identifiers are written out exactly as the sequencing file has them, including
the `t_` prefixes and `_NC_<date>` suffixes that re-sequenced samples carry,
because they have to match the FASTA headers. Those decorations are stripped
only to find the matching epidemiology row. When a specimen has several
sequences, the workflow keeps the one with the most unambiguous bases and
records the choice in `results/replicates.tsv`.

Read `conversion_report.txt`. It lists samples dropped for having no
epidemiology row and therefore no collection date, samples sequenced more than
once, and travel countries missing from the colour ordering file.

To write the table by hand instead, copy `defaults/metadata_template.tsv` and
fill one row per sample.

| Column | Required | Notes |
|---|---|---|
| `sample_id` | yes | must match a FASTA header, and must be unique |
| `serotype` | yes | `denv1` through `denv4`; `DENV1` and `1` are accepted and normalized |
| `nextclade_clade` | yes | the v-gen-lab lineage from Daytona, e.g. `2II_F.1.1.2` |
| `collection_date` | yes | `YYYY-MM-DD`, or `YYYY-MM-XX` / `YYYY-XX-XX` when partial |
| `location` | no | county |
| `case_origin` | no | `local`, `travel-associated`, or `undetermined`; set automatically when `travel_country` is filled |
| `travel_country` | no | where infection likely occurred; drives country_exposure |
| `host` | no | defaults to `Homo sapiens`; set to `Aedes aegypti` for vector pools |
| `strain` | no | overrides the derived Auspice display name |
| `country`, `region`, `division` | no | default to `USA`, `North America`, `Florida` |
| `authors`, `institution`, `notes` | no | default to the values in `defaults/config.yaml` |

Everything else is filled in automatically. `serotype_genbank`, `is_lab_host`,
`host_genus`, `host_type`, `length`, `data_source`, the three lineage levels, and
the two exposure columns are all derived, using the same values and spellings as
the public metadata so that colorings do not split into duplicate categories.

### Imported cases

Setting `travel_country` also sets `case_origin` to `travel-associated`, since
recording a travel country is what that means. The reverse is not inferred: a
blank `case_origin` with no travel history earns a warning rather than a default
of `local`, because an untravelled case and an uninvestigated one look identical
in the data and asserting local transmission you never established is the worse
error. Vector samples get a warning if `travel_country` is set, since mosquitoes
are collected where they are found.

Fill in `travel_country` and the workflow sets `country_exposure` to it, and
`region_exposure` to whichever region that country sits in according to
`phylogenetic/defaults/color_orderings.tsv`. Leave it blank and both fall back to
the collection country and region.

This matters because `augur traits` reconstructs ancestral geography from the
exposure columns, not from `country`. A case acquired in Cuba but reported in
Florida would otherwise make Florida look like the source of that lineage. The
map still places the sample in the USA, so it continues to read as Florida
surveillance; only the reconstruction uses the exposure value. This is the same
pattern the Nextstrain ncov builds used.

## Sample identifiers

Use whatever your lab already uses. `TVU26000019`, `JVV25001903`, and
`MosquitoPool_K26-10948` all work as-is.

The identifier lands in the `accession` column, which is what the phylogenetic
workflow keys tips on, and it must equal the FASTA header. Three properties
matter, and the validator checks all three:

- Only letters, digits, and `_ . -`. Spaces, slashes, colons, commas, and pipes
  break Newick trees and Auspice node names. This is `validate.id_regex` in
  `defaults/config.yaml`.
- Unique within a run, and stable across runs, so a sample keeps the same tip
  every time you rebuild.
- Not already a public GenBank accession. On a collision `augur merge` silently
  keeps one of the two sequences and discards the other, so the validator
  compares your identifiers against `../ingest/results/metadata_all.tsv` when
  that file exists.

The Auspice display name is derived for you as
`DENV2/USA/TVU26000019/2026`, matching the pattern the public GenBank records
use, so local and public tips read the same way in the tree. Supply a `strain`
column to override it.

## Samples already in GenBank

A genome the lab also deposited in GenBank comes back through `ingest` under a
GenBank accession, which the collision check above cannot see. The workflow finds
these copies by sequence. Each local genome is aligned with the GenBank records
from Florida of the same serotype. A pair counts as one genome when it agrees at
every position where both have a called base, and those positions cover at least
90% of the shorter sequence. When a genome matches more than one record, the one
with the same collection date is taken.

A linked sample takes the GenBank accession and drops its own sequence, so the
phylogenetic workflow shows it once: GenBank sequence and accession, with the
local date, case origin, travel country, county and `data_source`.
`results/genbank_copies.tsv` lists every pair. A pair marked `ambiguous` could
not be resolved and stays in the build twice, so decide those rows before
building. `bphl_named_unmatched` rows are GenBank records with a BPHL name and no
local genome, which is expected for samples absent from the local metadata.

## Running it

```sh
nextstrain build ingest
nextstrain build local
```

Run both from the top level of the repository. `local` reads `ingest`'s results
to find GenBank copies, so `ingest` has to finish first. Core count and Snakemake
flags come from [`profiles/default/config.yaml`](profiles/default/config.yaml),
so no `--cores` is needed.

Outputs land in `local/results/`:

- `metadata_{all,denv1..denv4}.tsv`
- `sequences_{all,denv1..denv4}.fasta`
- `include_{all,denv1..denv4}.txt`
- `validation_report.txt`
- `replicates.tsv`
- `genbank_copies.tsv`

All five serotypes always get a file, even when a run has no samples for one of
them, because the phylogenetic workflow expands its input paths over every
serotype on every run.

Read `validation_report.txt` before releasing a build. It lists per-sample
length, unambiguous base count, and ambiguity fraction, along with every error
and warning.

When validation fails, Snakemake deletes the outputs of the failed rule, so
`results/validation_report.txt` will not be there. The same errors and warnings
are in `logs/validate_local_metadata.txt`, which Snakemake leaves alone. Add
`--show-failed-logs` to have them printed to the terminal as well.

## What the validator rejects

The run stops, and nothing is written, on any of these:

- a required column missing, or an unrecognized column present
- an empty, duplicated, or malformed `sample_id`
- a serotype outside `denv1` to `denv4`
- an unparseable, future, or pre-1950 collection date
- `country` or `region` resolving to empty or `?`, which `augur filter` silently
  drops later
- a sample in the metadata with no FASTA record
- a duplicate FASTA header
- a character outside the IUPAC nucleotide alphabet

Warnings do not stop the run unless you set `validate.strict: true`. The ones
worth reading are short-sequence and high-ambiguity warnings, because those
samples enter the genome tree only because `include.txt` forces them past
`augur filter --min-length 5000`.

## On the lineage call

The `nextclade_clade` column comes straight from `Daytona_dengue`, which already
ran the `community/v-gen-lab/dengue` datasets when the sample was sequenced. This
workflow does not recompute it. It only splits the value into the three levels
Auspice colors by, using the same regex `ingest` uses:

```
2II_F.1.1.2  ->  genotype 2II,  major_lineage 2II_F,  minor_lineage 2II_F.1.1.2
```

The assumption is that Daytona and `ingest` called their lineages against the
same v-gen-lab dataset version. If v-gen-lab publishes a revision between your
sequencing run and your last `ingest`, your samples and the GenBank samples are
labelled by two different rulebooks and nothing here will warn you. Re-running
`ingest` around the same time you build keeps them aligned.

A sample whose `nextclade_clade` is `unclassified` or empty keeps its place in
the tree but shows blank under the lineage colorings, and the validator says so.

## A note on disclosure

The exported Auspice JSONs carry collection dates at day resolution and
county-level `location` values for local samples. Whether that is releasable
outside your organization is a policy question, not a technical one. Coarsen
`collection_date` to the month and drop `location` in the input file if the
answer is no.
