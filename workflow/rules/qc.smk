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


rule calculate_depth:
    input:
        fastq=rules.aggregate_pools.output.fastq,
        ref=infer_reference,
    output:
        depth=RESULTS / "depth/{experiment}.depth.tsv",
    log:
        LOGS / "calculate_depth/{experiment}.log",
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
