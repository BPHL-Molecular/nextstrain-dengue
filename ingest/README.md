# nextstrain-dengue/ingest

Fetches dengue sequences and metadata from NCBI GenBank, curates them, and
assigns dengue lineages with the `community/v-gen-lab/dengue` Nextclade
datasets. Output is the public backbone the phylogenetic build places local
sequences against.

## Software requirements

A Nextstrain runtime. Follow the [standard installation instructions](https://docs.nextstrain.org/en/latest/install.html),
then confirm with `nextstrain check-setup`.

## Usage

Run from the top level of the repository:

```sh
nextstrain build ingest
```

Fetch the raw NDJSON without running the rest:

```sh
nextstrain build ingest data/ncbi.ndjson
```

Core count and Snakemake flags come from [`profiles/default/config.yaml`](profiles/default/config.yaml),
so neither invocation needs `--cores`.

See the [HiPerGator section](../README.md#running-on-uf-hipergator) of the root
README for running this under `srun`, including what to do when compute nodes
cannot reach NCBI.

This produces 10 files (within the `ingest` directory):

A pair of files with all the dengue sequences:

- `ingest/results/metadata_all.tsv`
- `ingest/results/sequences_all.fasta`

A pair of files for each dengue serotype (denv1 - denv4)

- `ingest/results/metadata_denv1.tsv`
- `ingest/results/sequences_denv1.fasta`
- `ingest/results/metadata_denv2.tsv`
- `ingest/results/sequences_denv2.fasta`
- `ingest/results/metadata_denv3.tsv`
- `ingest/results/sequences_denv3.fasta`
- `ingest/results/metadata_denv4.tsv`
- `ingest/results/sequences_denv4.fasta`


## Configuration

Configuration takes place in `defaults/config.yaml` by default.

## Input data

### GenBank data

GenBank sequences and metadata are fetched via [NCBI datasets](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).
