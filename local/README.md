# local

Prepares privately sequenced dengue consensus genomes so the phylogenetic
workflow can place them on a tree alongside public GenBank data.

This workflow validates your metadata table, derives the dengue lineage levels
from the Nextclade call `Daytona_dengue` already made, and writes one metadata
file, one sequence file, and one include list per serotype.

## What you provide

Two things: the `Daytona_dengue` runs to take sequences from, and one metadata
table. Keep both under `local/input/`, which is gitignored so nothing private can
be committed by accident.

### The Daytona_dengue runs

List each run's output folder in a small YAML file, for example
`local/input/runs.yaml`:

```yaml
daytona_runs:
  - /blue/bphl-florida/share/daytona_output_2026
  - /blue/bphl-florida/share/daytona_output_2026_repeats
sample_metadata: input/metadata.txt
```

Each folder has to hold `summary_report.txt` and `assemblies_qc_pass/pass/` and
`assemblies_qc_pass/review/`, one consensus FASTA per sample. `Daytona_dengue`
writes every `PASS` assembly to `pass/` but only the `REVIEW` assemblies covering
at least 80% of the genome to `review/`, so the folders decide what is usable.
The workflow takes `sample_id`, `serotype`, `nextclade_clade` and `vadr_flag`
from the summary report and looks for a `PASS` sample only in `pass/` and a
`REVIEW` sample only in `review/`, by file name or by the identifier in the FASTA
header. A `REVIEW` sample absent from `review/` is left out as below the coverage
cutoff. Set `vadr_flags` to leave out `REVIEW` altogether. Samples with no serotype digit (`unclassified`, `NA`) are left out,
because there is no v-gen-lab dataset to place them against.

When the same `sample_id` appears in more than one run, a `PASS` copy is used
over a `REVIEW` one; between equal flags, the run listed last wins.

### `local/input/metadata.txt`

Tab-delimited, one row per sample, with these five columns:

| Column | Notes |
|---|---|
| `sample_id` | the laboratory identifier, without run tags |
| `collection_date` | `YYYY-MM-DD`, or `YYYY-MM-XX` / `YYYY-XX-XX` when partial |
| `location` | county, spelled as in `phylogenetic/defaults/lat_longs.tsv` (e.g. `Dade`) |
| `case_origin` | `local`, `travel-associated` or `undetermined` |
| `travel_country` | where infection likely occurred, for imported cases; drives `country_exposure` |

The table may list samples that were never sequenced or did not pass VADR; they
are ignored. A sequenced sample needs a row here, because the collection date
comes from it.

Run identifiers can carry tags the metadata does not: a `t_` prefix, a run suffix
(`_NC_<date>`, `_RJ_<date>`, `-repeat`, `_repeat`, `-repeat2`, `-NextSeq`,
`_test`) or a trailing `K` or `k`. They are stripped only to find the metadata
row; the tip keeps the full run identifier. When a specimen has several runs, the
workflow keeps a `PASS` run over any other, then the run with the most
unambiguous bases, and records the choice in `results/replicates.tsv`.

`travel_country` passes through `defaults/country_synonyms.tsv` to match the
spellings in `phylogenetic/defaults/color_orderings.tsv`. Mosquito pools are
recognized from `sample_id` (`vector_pattern`, default any identifier containing
"mosquito") and get `host` `Aedes aegypti` (`vector_host`), which the workflow
resolves to `Aedes` and `Mosquito` through ingest's host map.

### `results/input_report.txt`

Read it after every run. It lists how many samples each flag had, how many were
written, and every sample left out with the reason: no metadata row, a `PASS`
sample with no FASTA in `pass/`, a `REVIEW` sample below the coverage cutoff,
more than one matching FASTA, no serotype, sequenced in several runs, and
metadata rows that were not used.

The rest of each record is filled in automatically. `serotype_genbank`,
`is_lab_host`, `host_genus`, `host_type`, `length`, `data_source`, the three
lineage levels, and the two exposure columns are all derived, using the same
values and spellings as the public metadata so that colorings do not split into
duplicate categories. `country`, `region` and `division` default to `USA`,
`North America` and `Florida`.

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
nextstrain build local --configfile input/runs.yaml
```

Run both from the top level of the repository. `local` reads `ingest`'s results
to find GenBank copies, so `ingest` has to finish first. Core count and Snakemake
flags come from [`profiles/default/config.yaml`](profiles/default/config.yaml),
so no `--cores` is needed.

Outputs land in `local/results/`:

- `input_report.txt`
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
