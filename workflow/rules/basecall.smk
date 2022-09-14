rule basecall:
    input:
        fast5_dir=get_fast5_dir,
    output:
        summary=(
            RESULTS
            / f"basecalls/guppy_v{GUPPY_VERSION}/{{run}}/sequencing_summary.txt"
        ),
        save_path=directory(RESULTS / f"basecalls/guppy_v{GUPPY_VERSION}/{{run}}/"),
    log:
        LOGS / "basecall/{run}.log",
    threads: 2
    params:
        opt=" ".join(
            [
                f"-c {config['basecall_config']}",
                "--recursive",
                "--calib_detect",
                "--device cuda:all:100%",
            ]
        ),
        barcode_kits=get_barcode_kits,
    resources:
        mem_mb=6 * GB,
        time="1d",
        partition="gpgpu",
        slurm="gres=gpu:2 qos=gpgpumdhs",
    container:
        CONTAINERS["guppy"]
    shell:
        "guppy_basecaller {params.opt} {params.barcode_kits} -i {input.fast5_dir} -s {output.save_path} &> {log}"


rule merge_fastq:
    input:
        fastq_dir=rules.basecall.output.save_path,
    output:
        fastq=RESULTS / f"demux/guppy_v{GUPPY_VERSION}/{{run}}/{{sample}}.fq.gz",
    log:
        LOGS / "merge_fastq/{run}/{sample}.log",
    resources:
        mem_mb=GB,
        time="30m",
    params:
        rgx=r".*\.f(ast)?q(\.gz)?$",
        barcode_dir=get_barcode_dir,
    container:
        CONTAINERS["rs_utils"]
    resources:
        time="30m",
    shell:
        "fd -X cat \; {params.rgx:q} {input.fastq_dir}/pass/{params.barcode_dir} > {output.fastq} 2> {log}"
