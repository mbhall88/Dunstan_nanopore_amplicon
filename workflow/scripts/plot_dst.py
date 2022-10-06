import sys
from enum import Enum
from itertools import product

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


class Prediction(Enum):
    Resistant = "R"
    Susceptible = "S"
    MinorResistance = "r"
    Unknown = "U"
    Failed = "F"

    def __str__(self) -> str:
        return self.value


class Classification(Enum):
    TruePositive = "TP"
    FalsePositive = "FP"
    TrueNegative = "TN"
    FalseNegative = "FN"

    def __str__(self) -> str:
        return self.value


class Classifier:
    def __init__(
        self,
        minor_is_susceptible: bool = False,
        unknown_is_resistant: bool = False,
        failed_is_resistant: bool = False,
    ):
        self.minor_is_susceptible = minor_is_susceptible
        self.unknown_is_resistant = unknown_is_resistant
        self.failed_is_resistant = failed_is_resistant
        self.susceptible = {Prediction.Susceptible}
        self.resistant = {Prediction.Resistant}
        if self.minor_is_susceptible:
            self.susceptible.add(Prediction.MinorResistance)
        else:
            self.resistant.add(Prediction.MinorResistance)

        if self.unknown_is_resistant:
            self.resistant.add(Prediction.Unknown)
        else:
            self.susceptible.add(Prediction.Unknown)

        if self.failed_is_resistant:
            self.resistant.add(Prediction.Failed)
        else:
            self.susceptible.add(Prediction.Failed)

    def from_predictions(
        self, y_true: Prediction, y_pred: Prediction
    ) -> Classification:
        if y_true in self.susceptible:
            expected_susceptible = True
        elif y_true in self.resistant:
            expected_susceptible = False
        else:
            raise NotImplementedError(f"Don't know how to classify {y_true} calls yet")

        if y_pred in self.susceptible:
            called_susceptible = True
        elif y_pred in self.resistant:
            called_susceptible = False
        else:
            raise NotImplementedError(f"Don't know how to classify {y_pred} calls yet")

        if expected_susceptible and called_susceptible:
            return Classification.TrueNegative
        elif expected_susceptible and not called_susceptible:
            return Classification.FalsePositive
        elif not expected_susceptible and not called_susceptible:
            return Classification.TruePositive
        else:
            return Classification.FalseNegative


def main():
    phenotypes = pd.read_csv(snakemake.input.phenotypes)
    phenotypes = phenotypes.loc[~phenotypes["experiment"].isna()]
    phenotypes.set_index("experiment", drop=False, inplace=True, verify_integrity=True)

    df = pd.read_csv(snakemake.input.summary)
    df["drug"] = df["drug"].str.lower()
    experiments = set(df["experiment"])

    # this is redundant info from tbprofiler as it also lists all the drugs in the class
    skip_drugs = {"fluoroquinolones", "aminoglycosides", "all"}
    df.query("drug not in @skip_drugs and prediction != 'N'", inplace=True)

    drugs = set(df["drug"])

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

    df.set_index(["experiment", "drug"], inplace=True)
    # tbprofiler doesn't produce S calls
    for drug, exp in product(drug_order, experiments):
        ix = (exp, drug)
        if ix not in df.index:
            df.at[ix, "prediction"] = "S"

    df.reset_index(inplace=True, drop=False)
    preds = list(set(df["prediction"]))

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

    df = (
        df.sort_values(by="experiment")
        .sort_values(by="drug", key=sort_drugs)
        .reset_index(drop=True)
    )

    classifier = Classifier()
    disc = []
    for _, row in df.iterrows():
        exp = row["experiment"]
        drug = row["drug"]
        pred = Prediction(row["prediction"])
        try:
            ph = phenotypes.at[exp, drug]
        except KeyError:
            disc.append(False)
            continue
        if pd.isna(ph):
            disc.append(False)
            continue
        if len(ph) > 1:
            ph = "R"
        truth = Prediction(ph)
        clf = classifier.from_predictions(truth, pred)
        if clf in (Classification.FalseNegative, Classification.FalsePositive):
            disc.append(True)
        else:
            disc.append(False)

    df["discrepant"] = disc

    fig, ax = plt.subplots(figsize=FIGSIZE, dpi=DPI)
    sns.scatterplot(
        data=df,
        x="drug",
        y="experiment",
        hue="prediction",
        hue_order=sorted(preds, reverse=True),
        ax=ax,
        s=SIZE,
        edgecolor="black",
        palette=PALETTE,
        style="discrepant",
    )

    ax.set(xlabel="", ylabel="")
    xticks = [long2short[d] for d in drug_order]
    ax.set_xticklabels(xticks, fontsize=FS)
    ax.xaxis.tick_top()
    handles, labels = ax.get_legend_handles_labels()
    labels = [l.capitalize() for l in labels]
    ax.legend(
        handles=handles,
        labels=labels,
        bbox_to_anchor=(0.5, -0.05),
        ncol=2,
        prop=dict(size=FS),
        frameon=False,
        loc="lower center",
        markerscale=3,
    )
    _ = ax.tick_params(axis="y", which="major", labelsize=FS)

    plt.tight_layout()
    fig.savefig(snakemake.output.plot)


main()
