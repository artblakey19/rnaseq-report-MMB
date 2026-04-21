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
        gsea_h = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_hallmark.csv", contrast=CONTRAST_IDS),
        gsea_reactome = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_reactome.csv", contrast=CONTRAST_IDS),
        gsea_wikipathways = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_wikipathways.csv", contrast=CONTRAST_IDS),
        gsea_pid = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_pid.csv", contrast=CONTRAST_IDS),
        gsea_biocarta = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_biocarta.csv", contrast=CONTRAST_IDS),
        gsea_cgp = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_cgp.csv", contrast=CONTRAST_IDS),
        gsea_onco = expand(RESULTS / "enrichment" / "{contrast}" / "gsea_oncogenic.csv", contrast=CONTRAST_IDS),
        ora_gobp = expand(RESULTS / "enrichment" / "{contrast}" / "ora_gobp.csv", contrast=CONTRAST_IDS),
        ora_kegg = expand(RESULTS / "enrichment" / "{contrast}" / "ora_kegg.csv", contrast=CONTRAST_IDS),
        ora_reactome = expand(RESULTS / "enrichment" / "{contrast}" / "ora_reactome.csv", contrast=CONTRAST_IDS),
        ora_hallmark = expand(RESULTS / "enrichment" / "{contrast}" / "ora_hallmark.csv", contrast=CONTRAST_IDS),
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
