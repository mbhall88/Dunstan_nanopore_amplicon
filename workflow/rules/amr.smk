rule tbprofiler_predict:
    input:
        fastq=rules.aggregate_pools.output.fastq,
    output:
        report=RESULTS / "amr_predictions/tbprofiler/results/{experiment}.results.json",
    log:
        LOGS / "tbprofiler_predict/{experiment}.log",
    shadow:
        "shallow"
    resources:
        time="12h",
        mem_mb=lambda wildcards, attempt: attempt * int(12 * GB),
    threads: 4
    container:
        CONTAINERS["tbprofiler"]
    params:
        opts="--no_trim -p {experiment} --platform nanopore",
        outdir=lambda wildcards, output: Path(output.report).parent.parent,
    shell:
        "tb-profiler profile {params.opts} -1 {input.fastq} -t {threads} -d {params.outdir} &> {log}"


rule tbprofiler_collate:
    input:
        reports=expand(
            RESULTS / "amr_predictions/tbprofiler/results/{experiment}.results.json",
            experiment=EXPERIMENTS,
        ),
    output:
        outdir=directory(RESULTS / "reports/tbprofiler"),
    log:
        LOGS / "tbprofiler_collate.log",
    resources:
        time="30m",
    container:
        CONTAINERS["tbprofiler"]
    params:
        opts="--full",
        indir=lambda wildcards, input: Path(input.reports[0]).parent,
    shell:
        "tb-profiler collate {params.opts} -d {output.outdir} -p {output.outdir} &> {log}"
