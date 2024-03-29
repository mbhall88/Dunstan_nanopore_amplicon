import sys

sys.stderr = open(snakemake.log[0], "w")

import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use("ggplot")

samplesheet = snakemake.params.samplesheet
runs = samplesheet["runs"]
experiment = snakemake.wildcards.experiment
run, sample_id = experiment.split("_", maxsplit=1)
samples = runs[run]["samples"]

if "rpa" in sample_id:
    strategy = "RPA"
elif "pcr" in sample_id:
    strategy = "PCR"
else:
    m = re.search(r"Pool(?P<pool>\d+)", sample_id)
    if not m:
        _, sample = experiment.split("_", maxsplit=1)
        strategies = set()
        for s in samples:
            if s.startswith(sample) or s.endswith(sample):
                strategies.add(samples[s]["strategy"])
        if len(strategies) != 1:
            raise KeyError(f"Got more than one strategy for {experiment}")
        else:
            strategy = strategies.pop()
    else:
        pool = m.group("pool")
        strategy = samples[f"Pool{pool}"]["strategy"]

gene2pool = dict()
for pool, pool_genes in samplesheet["primer_pools"][strategy].items():
    if pool == "Pool16":
        continue
    for g in pool_genes:
        gene2pool[g] = pool

frames = []
for p in snakemake.input.depths:
    d = pd.read_csv(p, sep="\t", header=None, names=["gene", "pos", "depth"])
    d["time"] = Path(p).parts[-3]
    frames.append(d)


df = pd.concat(frames)
genes = sorted(set(df["gene"]))
print(df.groupby(["time", "gene"])["depth"].describe(), file=sys.stderr)

n_cells = len(genes)
if n_cells == 16:
    nrows = 4
    ncols = 4
elif n_cells == 24:
    nrows = 4
    ncols = 6
elif n_cells == 25:
    nrows = 5
    ncols = 5
elif n_cells == 12:
    nrows = 3
    ncols = 4
else:
    raise NotImplementedError(f"Don't know how many rows and cells for {n_cells} genes")

fig, axes = plt.subplots(
    nrows=nrows,
    ncols=ncols,
    figsize=(13, 13),
    dpi=300,
    sharex=True,
    sharey=True,
)

x = "time"
y = "depth"
CMAP = plt.get_cmap("Set2").colors
yticks = [
    int(i)
    for i in np.logspace(
        np.log10(max(1, df["depth"].min())), np.log10(df["depth"].max()), num=10
    )
]

is_pool16 = False

pool_nums = set()

for g, ax in zip(genes, axes.flatten()):
    if is_pool16:
        colour = CMAP[0]
    elif sample_id.startswith("Pool"):
        if "Pool16" in sample_id:
            colour = CMAP[0]
            is_pool16 = True
        else:
            pool = gene2pool[g]
            i = int(pool[-1]) - 1
            colour = CMAP[i]
            pool_nums.add(i + 1)
    else:
        for s, sample_info in samples.items():
            if sample_id in s:
                primers = sample_info["primers"]

                if primers == "Pool16":
                    colour = CMAP[0]
                    is_pool16 = True
                else:
                    pool = gene2pool[g]
                    i = int(pool[-1]) - 1
                    colour = CMAP[i]
                    pool_nums.add(i + 1)

                break

    data = df.query("gene==@g")
    sns.violinplot(data=data, x=x, y=y, ax=ax, color=colour, cut=0)
    ax.set(title=g, xlabel="Duration of sequencing")
    ax.set_yscale("log")

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    r = mpl.patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none", visible=False)


if not is_pool16:
    handles = []
    labels = []
    for i in sorted(pool_nums):
        c = CMAP[i - 1]
        pool = i
        label = f"Pool{pool}"
        line = mpl.patches.mlines.Line2D([], [], color=c, linewidth=5)
        labels.append(label)
        handles.append(line)
    fig.legend(
        handles,
        labels,
        bbox_to_anchor=(0.5, -0.03),
        loc="lower center",
        frameon=False,
        ncol=len(labels),
        prop=dict(size=12),
    )

plt.tight_layout()
fig.savefig(snakemake.output.plot, bbox_inches="tight")
