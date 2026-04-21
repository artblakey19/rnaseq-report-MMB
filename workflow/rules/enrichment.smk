rule gsea_hallmark:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_hallmark.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_hallmark.rds",
    params:
        collection = "H",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_hallmark.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule gsea_reactome:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_reactome.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_reactome.rds",
    params:
        collection = "C2:CP:REACTOME",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_reactome.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule gsea_wikipathways:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_wikipathways.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_wikipathways.rds",
    params:
        collection = "C2:CP:WIKIPATHWAYS",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_wikipathways.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule gsea_pid:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_pid.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_pid.rds",
    params:
        collection = "C2:CP:PID",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_pid.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule gsea_biocarta:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_biocarta.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_biocarta.rds",
    params:
        collection = "C2:CP:BIOCARTA",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_biocarta.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule gsea_cgp:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_cgp.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_cgp.rds",
    params:
        collection = "C2:CGP",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_cgp.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule gsea_oncogenic:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "gsea_oncogenic.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "gsea_oncogenic.rds",
    params:
        collection = "C6",
        ranking = config["enrichment"]["gsea"]["ranking"],
        min_size = config["enrichment"]["gsea"]["min_size"],
        max_size = config["enrichment"]["gsea"]["max_size"],
        seed = config["enrichment"]["gsea"]["seed"],
    log:
        "logs/enrichment/{contrast}_gsea_oncogenic.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/gsea.R"


rule ora_gobp:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "ora_gobp.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "ora_gobp.rds",
    params:
        database = "GO_BP",
        primary = config["de"]["primary"],
        secondary = config["de"]["secondary"],
        min_input_genes = config["enrichment"]["ora"]["min_input_genes"],
        max_input_genes = config["enrichment"]["ora"]["max_input_genes"],
    log:
        "logs/enrichment/{contrast}_ora_gobp.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/ora.R"


rule ora_kegg:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "ora_kegg.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "ora_kegg.rds",
    params:
        database = "KEGG",
        primary = config["de"]["primary"],
        secondary = config["de"]["secondary"],
        min_input_genes = config["enrichment"]["ora"]["min_input_genes"],
        max_input_genes = config["enrichment"]["ora"]["max_input_genes"],
    log:
        "logs/enrichment/{contrast}_ora_kegg.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/ora.R"


rule ora_reactome:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "ora_reactome.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "ora_reactome.rds",
    params:
        database = "Reactome",
        primary = config["de"]["primary"],
        secondary = config["de"]["secondary"],
        min_input_genes = config["enrichment"]["ora"]["min_input_genes"],
        max_input_genes = config["enrichment"]["ora"]["max_input_genes"],
    log:
        "logs/enrichment/{contrast}_ora_reactome.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/ora.R"


rule ora_hallmark:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        table = RESULTS / "enrichment" / "{contrast}" / "ora_hallmark.csv",
        rds = RESULTS / "enrichment" / "{contrast}" / "ora_hallmark.rds",
    params:
        database = "Hallmark",
        primary = config["de"]["primary"],
        secondary = config["de"]["secondary"],
        min_input_genes = config["enrichment"]["ora"]["min_input_genes"],
        max_input_genes = config["enrichment"]["ora"]["max_input_genes"],
    log:
        "logs/enrichment/{contrast}_ora_hallmark.log",
    conda:
        "../envs/r-enrichment.yaml"
    script:
        "../scripts/ora.R"


rule tf_activity:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
    output:
        scores = RESULTS / "tf_activity" / "{contrast}" / "tf_scores.tsv",
        top = RESULTS / "tf_activity" / "{contrast}" / "tf_top.tsv",
    params:
        method = config["tf_activity"]["method"],
        min_size = config["tf_activity"]["min_size"],
        padj_cutoff = config["tf_activity"]["padj_cutoff"],
        top_n = config["tf_activity"]["top_n"],
        split_complexes = config["tf_activity"]["split_complexes"],
    log:
        "logs/tf_activity/{contrast}.log",
    conda:
        "../envs/py-decoupler.yaml"
    script:
        "../scripts/tf_activity.py"


rule pathway_activity:
    input:
        de = RESULTS / "de" / "{contrast}" / "deseq2_results.csv",
        vst_matrix = RESULTS / "exploratory" / "vst_matrix.tsv",
    output:
        scores = RESULTS / "pathway_activity" / "{contrast}" / "progeny_scores.tsv",
    params:
        top_targets = config["pathway_activity"]["top_targets"],
    log:
        "logs/pathway_activity/{contrast}.log",
    conda:
        "../envs/py-decoupler.yaml"
    script:
        "../scripts/pathway_activity.py"
