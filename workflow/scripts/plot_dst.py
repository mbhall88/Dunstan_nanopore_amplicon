import sys

sys.stderr = open(snakemake.log[0], "w")

from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

plt.style.use("ggplot")

FS = snakemake.params.fontsize
DPI = snakemake.params.dpi
FIGSIZE = snakemake.params.figsize
PALETTE = snakemake.params.palette
SIZE = snakemake.params.marker_size


def main():
    df = pd.read_csv(snakemake.input.summary)
    df["drug"] = df["drug"].str.lower()

    # this is redundant info from tbprofiler as it also lists all the drugs in the class
    skip_drugs = {"fluoroquinolones", "aminoglycosides"}
    df.query("drug not in @skip_drugs", inplace=True)

    drugs = sorted(set(df["drug"]))
    preds = list(set(df["prediction"]))

    s = """AMK amikacin
    CAP capreomycin
    CFX ciprofloxacin
    DLM delamanid
    EMB ethambutol
    ETO ethionamide
    INH isoniazid
    KAN kanamycin
    LFX levofloxacin
    LZD linezolid
    MFX moxifloxacin
    OFX ofloxacin
    PZA pyrazinamide
    RIF rifampicin
    STM streptomycin"""
    long2short = dict()
    for line in s.splitlines():
        ab, d = line.split()
        long2short[d] = ab

    short2long = dict()
    for line in s.splitlines():
        ab, d = line.split()
        short2long[ab] = d

    first_line = [short2long[d] for d in ["INH", "RIF", "EMB", "PZA"]]
    fluoroquinolones = [short2long[d] for d in ["LFX", "MFX", "OFX", "CFX"]]  # group A
    macrolides = [
        short2long[d] for d in ["AMK", "CAP", "KAN", "STM"]
    ]  # group B - second-line injectables
    other = [short2long[d] for d in ["ETO", "LZD", "DLM"]]  # group C and D
    drug_order = [
        d for d in [*first_line, *fluoroquinolones, *macrolides, *other] if d in drugs
    ]

    def sort_drugs(a):
        xs = drug_order
        out = []
        c = Counter()
        for x in a:
            i = xs.index(x)
            d = xs[i]
            c[d] += 1
            out.append((i, c[d]))
        return out

    df = df.sort_values(by="drug", key=sort_drugs).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    sns.scatterplot(
        data=df,
        x="drug",
        y="run",
        hue="prediction",
        hue_order=sorted(preds, reverse=True),
        ax=ax,
        s=SIZE,
        edgecolor="black",
        palette=PALETTE,
    )

    ax.set(xlabel="", ylabel="")
    xticks = [long2short[d] for d in drug_order]
    ax.set_xticklabels(xticks, fontsize=FS)
    ax.xaxis.tick_top()
    handles, labels = ax.get_legend_handles_labels()
    pred_lookup = {"R": "Resistant", "r": "Minor resistance", "S": "Susceptible"}
    labels = [pred_lookup[l] for l in labels]
    ax.legend(
        handles=handles,
        labels=labels,
        bbox_to_anchor=(0.5, -0.05),
        ncol=len(preds),
        prop=dict(size=FS),
        frameon=False,
        loc="lower center",
        markerscale=3,
    )
    _ = ax.tick_params(axis="y", which="major", labelsize=FS)

    fig.savefig(snakemake.output.plot)


main()
