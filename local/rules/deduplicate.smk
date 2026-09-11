"""
This part of the workflow leaves one copy of every genome: replicate sequences
of a specimen are reduced to the best one, and samples already deposited in
GenBank are linked to that record instead of entering the build a second time.

Needs ingest/ to have finished, because the GenBank candidates come from its
results.

REQUIRED INPUTS:

    metadata  = data/metadata_lineages.tsv
    sequences = data/sequences_all.fasta

OUTPUTS:

    metadata  = data/metadata_linked.tsv
    sequences = data/sequences_linked.fasta
    reports   = results/replicates.tsv, results/genbank_copies.tsv
"""

LINKING_SEROTYPES = ["denv1", "denv2", "denv3", "denv4"]


rule resolve_replicates:
    input:
        metadata="data/metadata_lineages.tsv",
        sequences="data/sequences_all.fasta",
    output:
        metadata="data/metadata_dedup.tsv",
        sequences="data/sequences_dedup.fasta",
        report="results/replicates.tsv",
    log:
        "logs/resolve_replicates.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/resolve-replicates.py \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            --output-metadata {output.metadata:q} \
            --output-sequences {output.sequences:q} \
            --report {output.report:q}
        """


rule align_linking_candidates:
    input:
        metadata="data/metadata_dedup.tsv",
        sequences="data/sequences_dedup.fasta",
        public_metadata=config["genbank_copies"]["metadata"],
        public_sequences=config["genbank_copies"]["sequences"],
        reference=config["genbank_copies"]["reference"],
    output:
        alignment="data/linking/{serotype}/aligned.fasta",
    params:
        candidates="data/linking/{serotype}/candidates.fasta",
    wildcard_constraints:
        serotype="|".join(LINKING_SEROTYPES),
    threads: 4
    log:
        "logs/align_linking_candidates_{serotype}.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        augur filter \
            --metadata {input.metadata:q} \
            --metadata-id-columns accession \
            --sequences {input.sequences:q} \
            --query "serotype_genbank == '{wildcards.serotype}'" \
            --empty-output-reporting silent \
            --output-sequences {params.candidates:q}.local

        augur filter \
            --metadata {input.public_metadata:q} \
            --metadata-id-columns accession \
            --sequences {input.public_sequences:q} \
            --query "division == 'Florida'" \
            --query-columns division:str \
            --empty-output-reporting silent \
            --output-sequences {params.candidates:q}.public

        cat {params.candidates:q}.local {params.candidates:q}.public > {params.candidates:q}
        rm {params.candidates:q}.local {params.candidates:q}.public

        if [ -s {params.candidates:q} ]; then
            augur align \
                --sequences {params.candidates:q} \
                --reference-sequence {input.reference:q} \
                --output {output.alignment:q} \
                --fill-gaps \
                --remove-reference \
                --nthreads {threads}
        else
            : > {output.alignment:q}
        fi
        rm {params.candidates:q}
        """


rule link_genbank_copies:
    input:
        metadata="data/metadata_dedup.tsv",
        sequences="data/sequences_dedup.fasta",
        public_metadata=expand(config["genbank_copies"]["metadata"], serotype="all"),
        alignments=expand("data/linking/{serotype}/aligned.fasta", serotype=LINKING_SEROTYPES),
    output:
        metadata="data/metadata_linked.tsv",
        sequences="data/sequences_linked.fasta",
        report="results/genbank_copies.tsv",
    log:
        "logs/link_genbank_copies.txt",
    shell:
        r"""
        exec &> >(tee {log:q})

        python3 scripts/link-genbank-copies.py \
            --metadata {input.metadata:q} \
            --sequences {input.sequences:q} \
            --public-metadata {input.public_metadata:q} \
            --alignments {input.alignments:q} \
            --output-metadata {output.metadata:q} \
            --output-sequences {output.sequences:q} \
            --report {output.report:q}
        """
