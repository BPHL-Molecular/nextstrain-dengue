"""
This part of the workflow collects the phylogenetic tree and annotations to
export a Nextstrain dataset.
REQUIRED INPUTS:
    metadata        = results/{serotype}/metadata.tsv
    tree            = results/{serotype}/{gene}/tree.nwk
    branch_lengths  = results/{serotype}/{gene}/branch-lengths.json
    node_data       = results/{serotype}/{gene}/*.json
OUTPUTS:
    auspice_json      = auspice/dengue_{serotype}_{gene}.json
    frequencies_json  = auspice/dengue_{serotype}_{gene}_tip-frequencies.json
    There are optional sidecar JSON files that can be exported as part of the dataset.
    See Nextstrain's data format docs for more details on sidecar files:
    https://docs.nextstrain.org/page/reference/data-formats.html
This part of the workflow usually includes the following steps:
    - augur export v2
    - augur frequencies
See Augur's usage docs for these commands for more details.
"""

import json


rule colors:
    input:
        color_schemes = "defaults/color_schemes.tsv",
        color_orderings = "defaults/color_orderings.tsv",
        metadata = "results/{serotype}/metadata.tsv",
        manual_colors = "defaults/colors.tsv"
    output:
        colors = "results/{serotype}/colors.tsv"
    benchmark:
        "benchmarks/{serotype}/colors.txt"
    params:
        # Countries are grouped by region so each region gets its own full ramp
        # rather than a sliver of one shared interpolation. location has no
        # grouping to use, so it takes a single ramp.
        grouped_columns = ["country", "country_exposure"],
        flat_columns = ["location"],
    shell:
        """
        python3 scripts/assign-colors.py \
            --color-schemes {input.color_schemes} \
            --ordering {input.color_orderings} \
            --metadata {input.metadata} \
            --output {output.colors}.ramp

        # augur passes every row of the colours file through to the exported
        # scale without deduplicating, and Auspice then uses the first entry it
        # finds for a value. So the sections are concatenated in priority order,
        # manual pins first and the single shared ramp last, and the whole thing
        # is deduplicated on trait and value keeping the first row. That leaves
        # one colour per value and relies on no undocumented precedence.
        (
        cat {input.manual_colors}

        python3 scripts/assign-grouped-colors.py \
            --metadata {input.metadata} \
            --ordering {input.color_orderings} \
            --color-schemes {input.color_schemes} \
            --grouped-columns {params.grouped_columns} \
            --flat-columns {params.flat_columns}

        cat {output.colors}.ramp
        ) | awk -F'\t' '$0 == "" || !seen[$1, $2]++' > {output.colors}

        rm {output.colors}.ramp
        """


rule prepare_auspice_config:
    """Prepare the auspice config file for each serotypes"""
    output:
        auspice_config="results/defaults/{serotype}/{gene}/auspice_config.json",
    benchmark:
        "benchmarks/{serotype}/{gene}/prepare_auspice_config.txt"
    params:
        replace_clade_key=lambda wildcard: r"clade_membership" if wildcard.gene in ['genome'] else r"major_lineage",
        replace_clade_title=lambda wildcard: r"Serotype" if wildcard.serotype in ['all'] else r"Genotype (Nextclade)",
        # Auspice draws transmission lines by comparing a node's geographic value
        # with its parent's, so a geo resolution is only usable if augur traits
        # reconstructed it onto internal nodes. These must stay in step with
        # config.traits.traits_columns or the map silently draws nothing.
        geo_resolutions=lambda wildcard: (
            ["region", "region_exposure"]
            if wildcard.serotype == "all"
            else ["country", "region", "country_exposure", "region_exposure"]
        ),
    run:
        export_config = config.get("export", {})
        data = {
            "title": export_config.get("title", "Real-time tracking of dengue virus evolution"),
            "maintainers": export_config.get("maintainers", [
              {"name": "the Nextstrain team", "url": "https://nextstrain.org/team"}
            ]),
            "data_provenance": [
              {
                "name": "GenBank",
                "url": "https://www.ncbi.nlm.nih.gov/genbank/"
              }
            ],
            "build_url": export_config.get("build_url", "https://github.com/nextstrain/dengue"),
            "colorings": [
              {
                "key": "gt",
                "title": "Genotype",
                "type": "categorical"
              },
              {
                "key": "data_source",
                "title": "Data source",
                "type": "categorical"
              },
              {
                "key": "num_date",
                "title": "Date",
                "type": "continuous"
              },
              {
                "key": "country",
                "title": "Country",
                "type": "categorical"
              },
              {
                "key": "region",
                "title": "Region",
                "type": "categorical"
              },
              {
                "key": "serotype_genbank",
                "title": "Serotype (Genbank metadata)",
                "type": "categorical"
              },
              {
                "key": "genotype",
                "title": "Genotype (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "major_lineage",
                "title": "Major lineage (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "minor_lineage",
                "title": "Minor lineage (Nextclade)",
                "type": "categorical"
              },
              {
                "key": "host",
                "title": "Host Species",
                "type": "categorical"
              },
              {
                "key": "host_genus",
                "title": "Host Genus",
                "type": "categorical"
              },
              {
                "key": "host_type",
                "title": "Host Type",
                "type": "categorical"
              },
              {
                "key": "division",
                "title": "State or division",
                "type": "categorical"
              },
              {
                "key": "location",
                "title": "County",
                "type": "categorical"
              },
              {
                "key": "case_origin",
                "title": "Case origin",
                "type": "categorical"
              },
              {
                "key": "country_exposure",
                "title": "Country of exposure",
                "type": "categorical"
              },
              {
                "key": "region_exposure",
                "title": "Region of exposure",
                "type": "categorical"
              }
            ],
            "geo_resolutions": params.geo_resolutions,
            "display_defaults": {
              "map_triplicate": True,
              "color_by": params.replace_clade_key,
              "tip_label": "strain"
            },
            "filters": [
              "data_source",
              "case_origin",
              "location",
              "division",
              "country",
              "region",
              "author"
            ],
            "panels": [
              "tree",
              "map",
              "entropy",
              "frequencies"
            ],
            "metadata_columns": [
              "accession",
              "strain",
              "url",
              "data_source",
              "division",
              "location"
            ]
          }

        # During genome/dengue_all workflows, clade membership represents Serotype
        # While genome/dengue_denvX workflows, clade_membership represents the more detailed Genotype
        if params.replace_clade_key == 'clade_membership':
            if wildcards.gene in ['genome'] and wildcards.serotype in ['all']:
                clade_membership_title="Serotype (Nextstrain)"
            else:
                clade_membership_title="Genotype (Nextstrain)"

            data["colorings"].append({
                "key": "clade_membership",
                "title": clade_membership_title,
                "type": "categorical"
            })
        else:
            # During E/dengue_all workflows, default color by Serotype
            if wildcards.serotype in ['all']:
                data["display_defaults"]["color_by"]="serotype_genbank"

        with open(output.auspice_config, 'w') as fh:
            json.dump(data, fh, indent=2)


rule lat_longs:
    """
    augur keys coordinates by trait name as well as value, so country_exposure
    and region_exposure match nothing in its built-in table and augur export
    drops those geo resolutions entirely. This repeats the country and region
    rows under the exposure names. Every original row is kept, so it works
    whether --lat-longs supplements the built-in table or replaces it.

    defaults/lat_longs.tsv adds the countries augur's table does not carry.
    """
    input:
        extra = "defaults/lat_longs.tsv",
    output:
        lat_longs = "results/defaults/lat_longs.tsv",
    shell:
        r"""
        augur_lat_longs=$(python3 -c "import augur, pathlib; print(pathlib.Path(augur.__file__).parent / 'data' / 'lat_longs.tsv')")

        awk -F'\t' -v OFS='\t' '
            {{ print }}
            $1 == "country" {{ print "country_exposure", $2, $3, $4 }}
            $1 == "region"  {{ print "region_exposure",  $2, $3, $4 }}
        ' "$augur_lat_longs" {input.extra} > {output.lat_longs}
        """


rule export:
    """Exporting data files for auspice"""
    input:
        tree = "results/{serotype}/{gene}/tree.nwk",
        metadata = "results/{serotype}/metadata.tsv",
        branch_lengths = "results/{serotype}/{gene}/branch-lengths.json",
        traits = "results/{serotype}/{gene}/traits.json",
        clades = lambda wildcard: "results/{serotype}/{gene}/clades.json" if wildcard.gene in ['genome'] else [],
        nt_muts = "results/{serotype}/{gene}/nt-muts.json",
        aa_muts = "results/{serotype}/{gene}/aa-muts.json",
        description = config["export"]["description"],
        auspice_config = "results/defaults/{serotype}/{gene}/auspice_config.json",
        colors = "results/{serotype}/colors.tsv",
        lat_longs = "results/defaults/lat_longs.tsv",
    output:
        auspice_json = "auspice/dengue_{serotype}_{gene}.json"
    benchmark:
        "benchmarks/{serotype}/{gene}/export.txt"
    params:
        strain_id = config.get("strain_id_field", "strain"),
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --node-data {input.branch_lengths} {input.traits} {input.clades} {input.nt_muts} {input.aa_muts} \
            --colors {input.colors} \
            --lat-longs {input.lat_longs} \
            --description {input.description} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence-inline \
            --output {output.auspice_json}
        """

rule tip_frequencies:
    """
    Estimating KDE frequencies for tips
    """
    input:
        tree = "results/{serotype}/{gene}/tree.nwk",
        metadata = "results/{serotype}/metadata.tsv",
    output:
        tip_freq = "auspice/dengue_{serotype}_{gene}_tip-frequencies.json"
    benchmark:
        "benchmarks/{serotype}/{gene}/tip_frequencies.txt"
    params:
        strain_id = config["strain_id_field"],
        min_date = config["tip_frequencies"]["min_date"],
        max_date = config["tip_frequencies"]["max_date"],
        narrow_bandwidth = config["tip_frequencies"]["narrow_bandwidth"],
        wide_bandwidth = config["tip_frequencies"]["wide_bandwidth"]
    shell:
        r"""
        augur frequencies \
            --method kde \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.strain_id} \
            --min-date {params.min_date} \
            --max-date {params.max_date} \
            --narrow-bandwidth {params.narrow_bandwidth} \
            --wide-bandwidth {params.wide_bandwidth} \
            --output {output.tip_freq}
        """
