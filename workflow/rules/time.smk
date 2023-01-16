rule subset_by_time:
    input:
        fastq=rules.aggregate_pools.output.fastq,
    output:
        reads=RESULTS / "time/{time}/reads/{experiment}.{time}.fq.gz",
    log:
        LOGS / "subset_by_time/{time}/{experiment}.log",
    threads: 1
    resources:
        time="20m",
        mem_mb=lambda wildcards, attempt: attempt * int(1 * GB),
    container:
        CONTAINERS["ontime"]
    shell:
        "ontime --to {wildcards.time} -o {output.reads} {input.fastq} 2> {log}"


rule calculate_depth_for_time:
    input:
        fastq=rules.subset_by_time.output.reads,
        ref=infer_reference,
    output:
        depth=RESULTS / "time/{time}/depth/{experiment}.{time}.depth.tsv",
    log:
        LOGS / "calculate_depth_for_time/{time}/{experiment}.log",
    threads: 4
    resources:
        time="2h",
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
    conda:
        str(ENVS / "aln_tools.yaml")
    params:
        mm2_opts="-ax map-ont --secondary=no --sam-hit-only",
        samtools_opts="-aa",
    shell:
        """
        (minimap2 {params.mm2_opts} -t {threads} {input.ref} {input.fastq} \
            | samtools sort -@ {threads} \
            | samtools depth {params.samtools_opts} -o {output.depth} -) 2> {log}
        """

rule plot_depth_for_time:
    input:
        depths=expand(RESULTS / "time/{time}/depth/{{experiment}}.{time}.depth.tsv", time=TIMES),
    output:
        plot=report(
            PLOTS / "time/{experiment}.depth.png",
            category="Depth",
            subcategory="Time",
            labels={"experiment": "{experiment}", "time": "{time}"},
            ),
    log:
        LOGS / "plot_depth_for_time/{time}/{experiment}.log",
    conda:
        ENVS / "plot_depth.yaml"
    resources:
        time="15m",
    params:
        samplesheet=sample_data,
    script:
        str(SCRIPTS / "plot_time_depth.py")
