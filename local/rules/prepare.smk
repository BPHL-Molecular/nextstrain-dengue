"""
This part of the workflow builds the local metadata and sequences from
Daytona_dengue runs and the sample metadata.

REQUIRED INPUTS:

    runs     = config.daytona_runs, each with summary_report.txt and
               assemblies_qc_pass/{pass,review}/
    metadata = config.sample_metadata

OUTPUTS:

    metadata  = data/input_metadata.tsv
    sequences = data/input_sequences.fasta
    report    = results/input_report.txt
"""

import os


rule prepare_daytona_inputs:
    input:
        metadata=config["sample_metadata"],
        summary_reports=[os.path.join(run, "summary_report.txt") for run in config["daytona_runs"]],
        synonyms="defaults/country_synonyms.tsv",
    output:
        metadata="data/input_metadata.tsv",
        sequences="data/input_sequences.fasta",
        report="results/input_report.txt",
    params:
        runs=config["daytona_runs"],
        vadr_flags=config["vadr_flags"],
        vector_pattern=lambda wildcards: config["vector_pattern"],
        vector_host=config["vector_host"],
    log:
        "logs/prepare_daytona_inputs.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/prepare-daytona-inputs.py \
            --runs {params.runs:q} \
            --metadata {input.metadata:q} \
            --vadr-flags {params.vadr_flags:q} \
            --vector-pattern {params.vector_pattern:q} \
            --vector-host {params.vector_host:q} \
            --synonyms {input.synonyms:q} \
            --output-metadata {output.metadata:q} \
            --output-sequences {output.sequences:q} \
            --report {output.report:q}
        """
