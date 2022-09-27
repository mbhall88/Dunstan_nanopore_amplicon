rule plot_dst:
    input:
        summary=RESULTS / "reports/{tool}.csv",
        phenotypes=config["phenotypes"],
    output:
        plot=report(PLOTS / "{tool}.dst.png", category="DST prediction"),
    log:
        LOGS / "plot_dst/{tool}.log",
    resources:
        time="10m",
    params:
        fontsize=14,
        dpi=300,
        figsize=(13, 13),
        palette="Set2",
        marker_size=300,
    conda:
        ENVS / "plot_dst.yaml"
    script:
        str(SCRIPTS / "plot_dst.py")
