def get_fast5_dir(wildcards):
    return RUNS[wildcards.run]["run_dir"]


def get_barcode_kits(wildcards):
    kits = RUNS[wildcards.run]["barcode_kit"]
    s = " ".join(kits)
    return f'--barcode_kits "{s}"'

def get_barcode_dir(wildcards):
    barcode = RUNS[wildcards.run]["samples"][wildcards.sample]["barcode"]
    return barcode.replace("NB", "barcode")