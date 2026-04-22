rule drug_repositioning:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        hits = RESULTS / "drug_repositioning" / "{contrast}" / "l2s2_hits.tsv",
    params:
        top_up = config["drug_repositioning"]["top_up"],
        top_down = config["drug_repositioning"]["top_down"],
        service = config["drug_repositioning"]["service"],
    log:
        "logs/drug_repositioning/{contrast}.log",
    conda:
        "../envs/py-requests.yaml"
    script:
        "../scripts/l2s2_query.py"


rule render_report:
    input:
        qc_summary = RESULTS / "qc" / "qc_summary.tsv",
        pca = RESULTS / "exploratory" / "pca.rds",
        dds_vst = RESULTS / "exploratory" / "dds_vst.rds",
        de = expand(RESULTS / "de" / "{contrast}" / "deseq2_results.csv", contrast=CONTRAST_IDS),
        gsea = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_combined.csv", contrast=CONTRAST_IDS),
        ora = expand(RESULTS / "enrichment" / "{contrast}" / "ora_combined.csv", contrast=CONTRAST_IDS),
        tf = expand(RESULTS / "tf_activity" / "{contrast}" / "tf_scores.tsv", contrast=CONTRAST_IDS),
        progeny = expand(RESULTS / "pathway_activity" / "{contrast}" / "progeny_scores.tsv", contrast=CONTRAST_IDS),
        l2s2 = expand(RESULTS / "drug_repositioning" / "{contrast}" / "l2s2_hits.tsv", contrast=CONTRAST_IDS),
        template = "report/template.qmd",
    output:
        html = RESULTS / "report" / "report.html",
    params:
        config_path = "config/config.yaml",
        formats = config["report"]["formats"],
    log:
        "logs/report/render.log",
    conda:
        "../envs/quarto.yaml"
    script:
        "../scripts/render_report.R"
