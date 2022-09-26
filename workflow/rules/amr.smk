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
        mem_mb=lambda wildcards, attempt: attempt * int(30 * GB),
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
        """
        mkdir -p {output.outdir} 2> {log}
        tb-profiler collate {params.opts} -d {output.outdir} -p {output.outdir} &>> {log}
        """


rule mykrobe_predict:
    input:
        fastq=rules.aggregate_pools.output.fastq,
    output:
        report=RESULTS / "amr_predictions/mykrobe/results/{experiment}.results.json",
    log:
        LOGS / "mykrobe_predict/{experiment}.log",
    shadow:
        "shallow"
    resources:
        time="6h",
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
    threads: 2
    container:
        CONTAINERS["mykrobe"]
    params:
        opts=" ".join(
            [
                "--force",
                "-A",
                "-O json",
                "-D 0.20",
                "--species tb",
                "--sample {experiment}",
                "--ont",
            ]
        ),
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} \
            -i {input.fastq} -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """


rule combine_mykrobe_reports:
    input:
        reports=expand(str(RESULTS / "amr_predictions/mykrobe/results/{exp}.results.json"), exp=EXPERIMENTS),
    output:
        report=RESULTS / "amr_predictions/mykrobe/summary.csv",
    log:
        LOGS / "combine_mykrobe_reports.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_mykrobe_reports.py")
