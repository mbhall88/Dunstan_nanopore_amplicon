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
        report=RESULTS / "reports/tbprofiler.csv",
    log:
        LOGS / "tbprofiler_collate.log",
    resources:
        time="10m",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_tbprofiler_reports.py")


rule mykrobe_predict:
    input:
        fastq=rules.aggregate_pools.output.fastq,
        bed=infer_bed_file,
    output:
        report=RESULTS / "amr_predictions/mykrobe/results/{experiment}.results.json",
    log:
        LOGS / "mykrobe_predict/{experiment}.log",
    shadow:
        "shallow"
    resources:
        time="1d",
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
    threads: 8
    container:
        CONTAINERS["mykrobe"]
    params:
        opts=" ".join(
            [
                "--force",
                "-A",
                "-O json",
                "-D 0.10",
                "--species tb",
                "--sample {experiment}",
                "--ont",
            ]
        ),
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -T {input.bed}\
            -i {input.fastq} -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """


rule combine_mykrobe_reports:
    input:
        reports=expand(
            str(RESULTS / "amr_predictions/mykrobe/results/{exp}.results.json"),
            exp=EXPERIMENTS,
        ),
    output:
        report=RESULTS / "reports/mykrobe.csv",
    log:
        LOGS / "combine_mykrobe_reports.log",
    container:
        CONTAINERS["python"]
    script:
        str(SCRIPTS / "combine_mykrobe_reports.py")
