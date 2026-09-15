"""
This part of the workflow validates the metadata against the consensus FASTA and
normalises the sequence headers.

REQUIRED INPUTS:

    metadata  = data/input_metadata.tsv
    sequences = data/input_sequences.fasta

OUTPUTS:

    metadata  = data/metadata_validated.tsv
    sequences = data/sequences_all.fasta
    report    = results/validation_report.txt
"""

import json


rule validate_local_metadata:
    """
    Check the hand-filled TSV against the FASTA and emit a metadata table with
    the columns ingest produces. Fails the run on any structural problem.
    """
    input:
        metadata="data/input_metadata.tsv",
        sequences="data/input_sequences.fasta",
    output:
        metadata="data/metadata_validated.tsv",
        report="results/validation_report.txt",
    params:
        # These are lambdas because Snakemake expands wildcards in plain param
        # strings, and both the JSON blobs and the strain template contain
        # braces that are not wildcards.
        constants=lambda wildcards: json.dumps(config["constants"]),
        rules=lambda wildcards: json.dumps(config["validate"]),
        strain_template=lambda wildcards: config["strain_template"],
        strict="--strict" if config["validate"]["strict"] else "",
    log:
        "logs/validate_local_metadata.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/validate-local-metadata.py \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            --constants {params.constants:q} \
            --rules {params.rules:q} \
            --strain-template {params.strain_template:q} \
            --output-metadata {output.metadata:q} \
            --output-report {output.report:q} {params.strict}
        """


rule normalize_local_fasta:
    """
    Rewrite FASTA headers to the bare sample identifier, matching the bare
    accession headers ingest writes.
    """
    input:
        sequences="data/input_sequences.fasta",
        metadata="data/metadata_validated.tsv",
    output:
        sequences="data/sequences_all.fasta",
    log:
        "logs/normalize_local_fasta.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/normalize-local-fasta.py \
            --sequences {input.sequences:q} \
            --metadata {input.metadata:q} \
            --id-column accession \
            --output {output.sequences:q}
        """

