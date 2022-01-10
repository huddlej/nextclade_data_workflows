'''
This part of the workflow starts from files

  - pre-processed/sequences.fasta
  - pre-processed/metadata.tsv

and produces files

  - builds/{build_name}/sequences.fasta
  - builds/{build_name}/metadata.tsv

'''

localrules = sample_designated

build_dir = "builds"

rule prepare_build:
    input:
        sequences = build_dir + "/{build_name}/sequences.fasta",
        metadata = build_dir + "/{build_name}/metadata.tsv"

rule sample_designated:
    input:
        metadata = "pre-processed/open_pango_metadata.tsv",
    output:
        strains = "builds/pango_subsample.tsv"
    shell:
        """
        python scripts/sample_for_pango_build.py
        """

rule pango_sampling:
    input:
        sequences = "pre-processed/open_pango.fasta.xz",
        metadata = "pre-processed/open_pango_metadata.tsv",
        sample_strains = "builds/pango_subsample.tsv"
    output:
        sequences = build_dir + "/{build_name}/sequences_raw.fasta",
    log:
        "logs/subsample_{build_name}_{subsample}-pango.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude-all \
            --include {input.sample_strains} \
            --output {output.sequences} 2>&1 | tee {log}
        """

rule extract_metadata:
    input:
        strains = lambda w: [build_dir + f"/{w.build_name}/sample-{subsample}.txt"
                for subsample in config["builds"][w.build_name]["subsamples"]],
        metadata = "data/metadata.tsv"
    output:
        metadata = rules.prepare_build.input.metadata
    params:
        adjust = lambda w: config["builds"][w.build_name].get("metadata_adjustments",{}),
    benchmark:
        "benchmarks/extract_metadata_{build_name}.txt"
    run:
        import pandas as pd
        strains = set()
        for f in input.strains:
            with open(f) as fh:
                strains.update([x.strip() for x in fh if x[0]!='#'])

        d = pd.read_csv(input.metadata, index_col='strain', sep='\t').loc[list(strains)]
        if len(params.adjust):
            for adjustment  in params.adjust:
                ind = d.eval(adjustment["query"])
                d.loc[ind, adjustment['dst']] = d.loc[ind, adjustment['src']]

        d.to_csv(output.metadata, sep='\t')

rule exclude_outliers:
    input:
        sequences = "builds/{build_name}/sequences_raw.fasta",
        metadata = rules.prepare_build.input.metadata,
        exclude = "profiles/exclude.txt",
        sequence_index = "pre-processed/sequence_index.tsv",
    output:
        sampled_sequences = "builds/{build_name}/sequences.fasta",
        sampled_strains = "builds/{build_name}/subsample.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --sequence-index {input.sequence_index} \
            --exclude {input.exclude} \
            --output {output.sampled_sequences} \
            --output-strains {output.sampled_strains}
        """
