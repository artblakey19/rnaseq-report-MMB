log_con <- file(snakemake@log[[1]], open = "wt")
sink(log_con, type = "message")
sink(log_con, type = "output")

suppressPackageStartupMessages({
  library(dplyr)
  library(fgsea)
  library(msigdbr)
})

set.seed(42)

# --- Params ---------------------------------------------------------------
ranking <- snakemake@params[["ranking"]]
colls <- snakemake@params[["collections"]]

# --- Input / Output -------------------------------------------------------
de_path <- snakemake@input[["de"]]
table_out <- snakemake@output[["table"]]

de_res <- read.csv(de_path, header = TRUE, stringsAsFactors = FALSE)

# Calculate gene ranks
if (ranking == "stat") {
  ranks_df <- de_res %>%
    filter(!is.na(gene_id), !is.na(stat)) %>%
    transmute(gene_id, metric = stat) %>%
    arrange(desc(metric))
} else if (ranking == "signed_p") {
  ranks_df <- de_res %>%
    filter(!is.na(gene_id), !is.na(log2FoldChange), !is.na(pvalue)) %>%
    transmute(gene_id, metric = sign(log2FoldChange) * (-log10(pmax(pvalue, 1e-300)))) %>%
    arrange(desc(metric))
} else {
  stop("Unsupported ranking method: ", ranking)
}
# fgsea requires named vector
ranks <- setNames(ranks_df$metric, ranks_df$gene_id)

# --- MSigDB Data Prep & GSEA ----------------------------------------------
get_pathways <- function(coll_str) {
  parts <- strsplit(coll_str, ":", fixed = TRUE)[[1]]
  m <- msigdbr(
    species = "Homo sapiens",
    collection = parts[1],
    subcollection = if (length(parts) > 1) paste(parts[-1], collapse = ":") else NULL
  )
  split(m$ensembl_gene, m$gs_name)
}

run_gsea <- function(coll) {
  message("Running GSEA for collection: ", coll)
  res <- fgsea(get_pathways(coll), ranks, minSize = 15, maxSize = 500)
  if (!nrow(res)) {
    return(NULL)
  }
  res$collection <- coll
  res$leadingEdge <- vapply(res$leadingEdge, paste, character(1), collapse = ";")
  as.data.frame(res)
}

col_order <- c("collection", "pathway", "pval", "padj", "log2err", "ES", "NES", "size", "leadingEdge")
# Run GSEA for each collection
combined_df <- bind_rows(lapply(colls, run_gsea))
combined_df <- combined_df[, intersect(col_order, colnames(combined_df))]

dir.create(dirname(table_out), showWarnings = FALSE, recursive = TRUE)
write.csv(combined_df, table_out, row.names = FALSE)
