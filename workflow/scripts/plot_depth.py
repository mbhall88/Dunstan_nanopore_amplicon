import sys
sys.stderr = open(snakemake.log[0], "w")

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

plt.style.use("ggplot")

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
    data = df.query("gene==@g")
    sns.lineplot(data=data, x=x, y=y, ax=ax)
    ax.set(title=g, xlabel="Relative position")
    ax.set_yscale("log")

    ax.set_yticks(yticks)
    ax.set_yticklabels(yticks)
    r = mpl.patches.Rectangle((0, 0), 1, 1, fill=False, edgecolor="none", visible=False)
    ax.legend([r], [f"median depth = {int(median_depth[g])}"], frameon=False)

fig.savefig(snakemake.output.plot)
