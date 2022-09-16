rule qc_plot:
    input:
        fastq=rules.merge_fastq.output.fastq,
    output:
        outdir=directory(RESULTS / "QC/plots/{run}/{sample}"),
    log:
        LOGS / "qc_plot/{run}/{sample}.log",
    threads: 2
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(8 * GB),
        time="3h",
    params:
        opts=" ".join(
            [
                "-p {sample}",
                "--tsv_stats",
                "--loglength",
                "--readtype 1D",
                "--no-N50",
                "--title 'Run: {run} Sample: {sample}'",
                "--dpi 300",
            ]
        ),
    container:
        CONTAINERS["nanoplot"]
    shell:
        "NanoPlot -t {threads} {params.opts} --fastq_rich {input.fastq} -o {output.outdir} &> {log}"


rule aggregate_pools:
    input:
        fastqs=infer_fastqs_to_aggregate,
    output:
        fastq=RESULTS / "combined_reads/{experiment}.fq.gz",
    log:
        LOGS / "aggregate_pools/{experiment}.log",
    resources:
        time="20m",
    container:
        CONTAINERS["rs_utils"]
    shell:
        "cat {input.fastqs} > {output.fastq} 2> {log}"