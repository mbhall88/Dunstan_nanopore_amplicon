import re


def get_fast5_dir(wildcards):
    return RUNS[wildcards.run]["run_dir"]


def get_barcode_kits(wildcards):
    kits = RUNS[wildcards.run]["barcode_kit"]
    s = " ".join(kits)
    return f'--barcode_kits "{s}"'


def get_barcode_dir(wildcards):
    barcode = RUNS[wildcards.run]["samples"][wildcards.sample]["barcode"]
    return barcode.replace("NB", "barcode")


def infer_fastqs_to_aggregate(wildcards):
    exp = wildcards.experiment
    fastqs = []
    run, sample_id = exp.split("_", maxsplit=1)
    if sample_id.startswith("Pool"):
        m = re.search(r"Pool(?P<start>\d+)-?(?P<end>\d+)?", sample_id)
        if not m:
            raise ValueError(f"Got unknown experiment {exp}")
        start = int(m.group("start"))
        end = m.group("end")
        if end is None:
            pools = [start]
        else:
            pools = list(range(start, int(end) + 1))
    else:
        pools = [1, 2, 3]

    for p in pools:
        if sample_id.startswith("Pool"):
            s = f"Pool{p}"
        else:
            s = f"Pool{p}-{sample_id}"

        fastqs.append(RESULTS / f"demux/guppy_v{GUPPY_VERSION}/{run}/{s}.fq.gz")

    return fastqs
