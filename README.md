# nextstrain-dengue

Dengue virus phylogenetics for the Florida Bureau of Public Health Laboratories,
forked from [nextstrain/dengue](https://github.com/nextstrain/dengue).

This fork runs entirely on local infrastructure. It fetches public GenBank data
for global context, adds consensus genomes from the
[Daytona_dengue](https://github.com/BPHL-Molecular/Daytona_dengue) pipeline, and
writes Auspice JSON files that colleagues open by drag and drop at
[auspice.us](https://auspice.us). Nothing reads from or writes to AWS.

## Workflows

Run them in this order.

1. [`ingest/`](./ingest) downloads dengue sequences from NCBI GenBank, curates
   them, and assigns dengue lineages with the `community/v-gen-lab/dengue`
   Nextclade datasets. Re-run it when you want fresher public context.
2. [`local/`](./local) validates your metadata, derives the lineage levels from
   the Nextclade call Daytona already made, and writes per-serotype files. Re-run
   it whenever you have new samples.
3. [`phylogenetic/`](./phylogenetic) merges the two and builds the trees.

[`nextclade/`](./nextclade) builds the Nextstrain genotype datasets and is
inherited from upstream. It is not part of the Florida build.

## Output

The Florida build produces ten Auspice datasets, the full cross product of five
serotype groupings and two genes:

| | Whole genome | E gene |
|---|---|---|
| All serotypes together | `dengue_all_genome.json` | `dengue_all_E.json` |
| DENV1 | `dengue_denv1_genome.json` | `dengue_denv1_E.json` |
| DENV2 | `dengue_denv2_genome.json` | `dengue_denv2_E.json` |
| DENV3 | `dengue_denv3_genome.json` | `dengue_denv3_E.json` |
| DENV4 | `dengue_denv4_genome.json` | `dengue_denv4_E.json` |

Each comes with a `_tip-frequencies.json` sidecar. All twenty land in
`phylogenetic/auspice/`.

Every tree carries the dengue lineage system of Hill et al. 2024 at three levels
of detail, from the Nextclade datasets at
[community/v-gen-lab/dengue](https://github.com/nextstrain/nextclade_data/tree/master/data/community/v-gen-lab/dengue),
plus a "Data source" coloring that separates local sequences from public ones.

## Quick start

Everything runs through the [Nextstrain CLI](https://docs.nextstrain.org/projects/cli/),
which supplies the toolchain from its own managed runtime. Confirm you have one
with `nextstrain check-setup`, then, from the repo root:

Put the sample metadata at `local/input/metadata.txt` and list the
`Daytona_dengue` output folders in `local/input/runs.yaml` (see
[`local/README.md`](local/README.md)), then:

```sh
nextstrain build ingest
nextstrain build local --configfile input/runs.yaml
nextstrain build phylogenetic --configfile build-configs/florida/config.yaml --cores 8
```

`ingest` and `local` each carry a `profiles/default/`, so Snakemake picks up
their core count and flags on its own. `phylogenetic` has none, matching
upstream, so give it `--cores` explicitly.

Run them in this order: `local` reads `ingest`'s results to link samples already
in GenBank, and `phylogenetic` needs both.

## Running on UF HiPerGator

Conda is the only usable runtime here: docker is not installed, and singularity
fails because the administrators disabled overlay support, which
`--writable-tmpfs` requires. `nextstrain check-setup` will confirm that.

Get an allocation, then work inside it:

```sh
srun --account=bphl-umbrella --qos=bphl-umbrella \
     --cpus-per-task=8 --mem=32gb --time=08:00:00 --pty bash -i

conda activate nextstrain          # provides the `nextstrain` command
cd /blue/bphl-florida/$USER/nextstrain-dengue
nextstrain build ingest
```

Rough sizing per stage, to tune afterwards from `ingest/benchmarks/` and
`phylogenetic/benchmarks/`:

| Stage | Cores | Memory | Time |
|---|---|---|---|
| `ingest` | 8 | 32 GB | a few hours, mostly the NCBI download |
| `local` | 4 | 8 GB | minutes |
| `phylogenetic` | 32 | 128 GB | long |

Before the phylogenetic build, pin the linear algebra libraries to one thread
each:

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
nextstrain build phylogenetic --configfile build-configs/florida/config.yaml --cores 32
```

`augur align` and `augur tree` already take eight threads apiece from Snakemake.
Without those exports, the numpy underneath several concurrent `augur refine`
runs will oversubscribe the allocation and slow everything down.

If one allocation cannot cover all ten builds, run a serotype at a time by naming
its two targets:

```sh
nextstrain build phylogenetic --configfile build-configs/florida/config.yaml --cores 8 \
    auspice/dengue_denv1_genome.json auspice/dengue_denv1_E.json
```

The `all` build is the long pole, since it carries every serotype.

If your partition's compute nodes have no outbound network, fetch the
network-dependent targets on a login node first:

```sh
nextstrain build ingest --notemp \
    data/ncbi_dataset.zip \
    data/nextclade_data/v-gen-lab/denv{1,2,3,4}.zip
```

`--notemp` is not optional. The fetch rule marks `data/ncbi_dataset.zip` as a
temporary file, so without it Snakemake deletes the download as soon as the rules
consuming it finish.

## Private data

`local/input/`, `local/data/`, and `local/results/` are all gitignored, so
consensus sequences and patient metadata cannot be committed by accident. Keep
it that way.

## Documentation

- [Contributor documentation](./CONTRIBUTING.md)
