import sys

sys.stderr = open(snakemake.log[0], "w")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import re
import seaborn as sns

plt.style.use("ggplot")

samplesheet = snakemake.params.samplesheet
experiment = snakemake.wildcards.experiment
run, sample_id = experiment.split("_", maxsplit=1)
if "rpa" in sample_id:
    strategy = "RPA"
elif "pcr" in sample_id:
    strategy = "PCR"
else:
    samples = samplesheet["runs"][run]["samples"]
    m = re.search(r"Pool(?P<pool>\d+)", sample_id)
    if not m:
        raise ValueError(f"Can't infer pool from {experiment}")
    pool = m.group("pool")
    strategy = samples[f"Pool{pool}"]["strategy"]
gene2pool = dict()
for pool, pool_genes in samplesheet["primer_pools"][strategy].items():
    if pool == "Pool16":
        continue
    for g in pool_genes:
        gene2pool[g] = pool

df = pd.read_csv(
    snakemake.input.depth, sep="\t", header=None, names=["gene", "pos", "depth"]
)
df["relpos"] = df.groupby("gene")["pos"].transform(lambda xs: xs / xs.max())
df["reldepth"] = df.groupby("gene")["depth"].transform(lambda xs: xs / xs.median())
genes = sorted(set(df["gene"]))
median_depth = dict()
for g in genes:
    median_depth[g] = df.query("gene==@g")["depth"].median()

fig, axes = plt.subplots(
    nrows=4,
    ncols=4,
    figsize=(13, 13),
    dpi=300,
    sharex=True,
    sharey=True,
    constrained_layout=True,
)

x = "relpos"
y = "depth"
CMAP = plt.get_cmap("Set2").colors
yticks = [
    int(i)
    for i in np.logspace(
        np.log10(max(1, df["depth"].min())), np.log10(df["depth"].max()), num=10
    )
]

for g, ax in zip(genes, axes.flatten()):
    if "Pool16" in sample_id:
        colour = CMAP[0]
    else:
        pool = gene2pool[g]
        i = int(pool[-1]) - 1
        colour = CMAP[i]

    data = df.query("gene==@g")
    sns.lineplot(data=data, x=x, y=y, ax=ax, color=colour)
    ax.set(title=g, xlabel="Relative position")
    ax.set_yscale("log")

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    r = mpl.patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none", visible=False)
    ax.legend([r], [f"median depth = {int(median_depth[g])}"], frameon=False)

if "Pool16" not in sample_id:
    handles = []
    labels = []
    for i in range(3):
        c = CMAP[i]
        pool = i + 1
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

fig.savefig(snakemake.output.plot)
