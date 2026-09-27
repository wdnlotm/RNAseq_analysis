library(patchwork)
library(ggplot2)
library(DESeq2)   # for plotCounts()
library(msigdbr)

#' Plot normalized counts for a differentially expressed gene
#'
#' @param ii Integer index of the row in `genes` to plot.
#' @param genes A data frame with `ensembl_id` and `hgnc_symbol` columns.
#' @param dds A DESeqDataSet.
#' @param group Name of the `colData` variable to plot on the x-axis.
#'
#' @return A ggplot object.
#' @examples
#' plot_gene_counts(1, diff_genes, dds, group = "etf_12")
plot_gene_counts <- function(ii, genes, dds, group) {
  ens <- genes$ensembl_id[ii]
  genename <- genes$hgnc_symbol[ii]
  # fall back to the Ensembl ID if the symbol is missing/empty
  if (is.na(genename) || genename == "") {
    genename <- ens
  }
  
  d <- plotCounts(dds, gene = ens, intgroup = group, returnData = TRUE)
  
  ggplot(d, aes(x = .data[[group]], y = count)) +
    geom_point(
      position = position_jitter(width = 0.1, height = 0),
      aes(color = .data[[group]])
    ) +
    scale_y_log10(breaks = c(25, 100, 400)) +
    ggtitle(genename)
}



boxplot_gene_counts <- function(ii, genes, dds, group) {
  ens <- genes$ensembl_id[ii]
  genename <- genes$hgnc_symbol[ii]
  # fall back to the Ensembl ID if the symbol is missing/empty
  if (is.na(genename) || genename == "") {
    genename <- ens
  }
  
  d <- plotCounts(dds, gene = ens, intgroup = group, returnData = TRUE)
  
  ggplot(d, aes(x = .data[[group]], y = count)) +
    geom_point(position = position_jitter(w = 0.1, h = 0), aes(color = .data[[group]])) +
    scale_y_log10(breaks = c(25, 100, 400)) +
    ggtitle(genename[1])
  
}



#' Run GO enrichment on a set of differentially expressed genes
#'
#' @param diff_genes A data frame with an `hgnc_symbol` column.
#' @param onto GO ontology to test: "BP", "MF", "CC", or "ALL".
#' @param universe Character vector of Entrez IDs to use as the enrichment
#'   background, normally the genes tested by DESeq2. `NULL` falls back to
#'   clusterProfiler's default of all annotated genes in `orgdb`, which is
#'   usually too broad and biases results toward expressed-gene terms.
#' @param orgdb An OrgDb annotation object.
#' @param pvalue_cutoff Adjusted p-value cutoff passed to `enrichGO()`.
#' @param bar_count Maximum number of categories to show in the barplot.
#' @param out_dir Directory for the results CSV. Set to `NULL` to skip writing.
#' @param plot Whether to draw a barplot of the enriched terms.
#'
#' @return The `enrichResult` object, invisibly, or `NULL` if nothing was found.
run_go_enrichment <- function(diff_genes,
                              onto,
                              universe = NULL,
                              orgdb = org.Hs.eg.db::org.Hs.eg.db,
                              pvalue_cutoff = PVALUECUTOFF,
                              qvalue_cutoff = 0.2,
                              bar_count = ENRICH_BARCOUNT,
                              out_dir = ".",
                              plot = TRUE) {
  
  converted_genes <- clusterProfiler::bitr(
    geneID   = diff_genes$hgnc_symbol,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = orgdb
  )
  
  if (nrow(converted_genes) == 0) {
    message("No symbols could be mapped to Entrez IDs")
    return(invisible(NULL))
  }
  
  if (is.null(universe)) {
    warning("No `universe` supplied: enrichment is tested against the whole ",
            "OrgDb rather than the genes tested in this experiment.",
            call. = FALSE)
  } else {
    universe <- unique(universe)
    n_outside <- sum(!converted_genes$ENTREZID %in% universe)
    if (n_outside > 0) {
      stop(n_outside, " DE genes are absent from `universe`; the background ",
           "does not match the tested gene set.", call. = FALSE)
    }
  }
  
  go_results <- clusterProfiler::enrichGO(
    gene          = converted_genes$ENTREZID,
    universe      = universe,
    OrgDb         = orgdb,
    keyType       = "ENTREZID",
    ont           = onto,
    pvalueCutoff  = pvalue_cutoff,
    qvalueCutoff  = qvalue_cutoff,
    pAdjustMethod = "BH",
    readable      = TRUE
  )
  
  n_terms <- if (is.null(go_results)) 0L else nrow(as.data.frame(go_results))
  
  if (n_terms == 0) {
    message("No enrichGO result for ontology ", onto)
    return(invisible(go_results))
  }
  
  if (!is.null(out_dir)) {
    readr::write_csv(
      as.data.frame(go_results),
      file.path(out_dir, paste0("go_df_", onto, ".csv"))
    )
  }
  
  if (plot) {
    print(barplot(go_results, showCategory = min(n_terms, bar_count), font.size = 9))
  }
  
  invisible(go_results)
}


#' Run fgsea against an MSigDB collection and build summary plots
#'
#' @param category MSigDB collection, e.g. "H", "C2", "C5".
#' @param gene_ranks Named numeric vector of gene-level statistics, keyed by
#'   Ensembl gene ID and sorted or unsorted (fgsea handles ordering).
#' @param n_top Number of up- and down-regulated pathways for the summary table.
#' @param n_enrich Number of up- and down-regulated pathways to plot
#'   individually. Capped at `n_top`.
#' @param min_size,max_size Pathway size limits passed to `fgsea()`.
#' @param padj_cutoff Adjusted p-value cutoff for pathway selection. Use 1 to
#'   rank purely by p-value, as in the original.
#' @param gsea_param Weighting exponent, used for both the table plot and the
#'   individual enrichment curves so they are comparable.
#' @param seed Seed for fgsea's randomized sampling. `NULL` to skip.
#' @param bpparam BiocParallel backend for `fgsea()`. Defaults to `SerialParam()`
#'   because the Windows default (`SnowParam`, a PSOCK cluster) emits
#'   "'package:stats' may not be available when loading" warnings while
#'   serializing to workers, and the cluster startup cost exceeds the gain for
#'   collections of this size. Pass e.g. `BiocParallel::MulticoreParam()` on
#'   Linux/macOS if a collection is large enough to need it.
#'
#' @return A named list with `results`, `table_plot`, and `enrichment_plots`.
msigdb_function <- function(category,
                            gene_ranks,
                            n_top = 10,
                            n_enrich = 2,
                            min_size = 15,
                            max_size = 500,
                            padj_cutoff = 1,
                            gsea_param = 0.5,
                            seed = 4827,
                            bpparam = BiocParallel::SerialParam()) {
  
  n_enrich <- min(n_enrich, n_top)
  
  # msigdbr renamed `category` to `collection` in recent versions
  msig_args <- list(species = "Homo sapiens")
  if ("collection" %in% names(formals(msigdbr::msigdbr))) {
    msig_args$collection <- category
  } else {
    msig_args$category <- category
  }
  msig_human <- do.call(msigdbr::msigdbr, msig_args)
  human_pathway_list <- split(msig_human$ensembl_gene, msig_human$gs_name)
  
  if (!is.null(seed)) set.seed(seed)
  
  fgsea_results <- fgsea::fgsea(
    pathways = human_pathway_list,
    stats    = gene_ranks,
    minSize  = min_size,
    maxSize  = max_size,
    BPPARAM  = bpparam
  )
  
  # plain data.frame copy, so selection does not depend on data.table NSE
  res <- as.data.frame(fgsea_results)
  
  pick <- function(direction) {
    keep <- !is.na(res$NES) & !is.na(res$padj) & res$padj <= padj_cutoff &
      if (direction > 0) res$NES > 0 else res$NES < 0
    hits <- res[keep, , drop = FALSE]
    head(hits$pathway[order(hits$pval)], n_top)
  }
  
  up   <- pick(1)
  down <- pick(-1)
  
  if (length(up) == 0 && length(down) == 0) {
    message("No pathways passed the size and significance filters for ", category)
    return(invisible(list(results = fgsea_results,
                          table_plot = NULL,
                          enrichment_plots = NULL)))
  }
  
  top_pathways <- c(up, rev(down))
  
  table_plot <- fgsea::plotGseaTable(
    human_pathway_list[top_pathways],
    gene_ranks,
    fgsea_results,
    pathwayLabelStyle = list(size = 7),
    valueStyle        = list(size = 8),
    headerLabelStyle  = list(size = 8),
    gseaParam         = gsea_param
  )
  
  # strongest pathways in each direction, by name rather than by position
  enrich_names <- c(head(up, n_enrich), head(down, n_enrich))
  
  enrichment_plots <- lapply(enrich_names, function(nm) {
    nes <- res$NES[match(nm, res$pathway)]
    padj <- res$padj[match(nm, res$pathway)]
    fgsea::plotEnrichment(human_pathway_list[[nm]], gene_ranks,
                          gseaParam = gsea_param) +
      labs(title = paste(strwrap(gsub("_", " ", nm), 40), collapse = "\n"),
           subtitle = sprintf("%s | NES = %.2f | padj = %.3g",
                              if (nes > 0) "Up" else "Down", nes, padj)) +
      theme(plot.title = element_text(size = 10),
            plot.subtitle = element_text(size = 8))
  })
  
  list(
    results          = fgsea_results,
    table_plot       = table_plot,
    enrichment_plots = patchwork::wrap_plots(enrichment_plots, ncol = 2)
  )
}


plotDists <- function(vsd.obj, label_col = NULL,  fontsize = 8) {
  sampleDists <- dist(t(assay(vsd.obj)))
  sampleDistMatrix <- as.matrix(sampleDists)

  if (!is.null(label_col)) {
    stopifnot(label_col %in% colnames(colData(vsd.obj)))
    # labels <- paste(colnames(vsd.obj), as.character(colData(vsd.obj)[[label_col]]), sep = " | ")
    rownames(sampleDistMatrix) <- paste(as.character(colData(vsd.obj)[[label_col]]))
    colnames(sampleDistMatrix) <- paste(colnames(vsd.obj))
  }

  colors <- colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Blues")))(255)
  pheatmap::pheatmap(
    sampleDistMatrix,
    clustering_distance_rows = sampleDists,
    clustering_distance_cols = sampleDists,
    col = colors,
    fontsize = fontsize,
    fontsize_row = fontsize-2,
    fontsize_col = fontsize,
    # cellwidth  = 7,
    # cellheight = 7   # equal values give square cells
  )
}


variable_gene_heatmap <- function(vsd.obj,
                                  num_genes = 500,
                                  annotation,
                                  title = "",
                                  label_cutoff = 60) {
  stopifnot(
    !missing(annotation),
    c("ensembl_gene_id", "hgnc_symbol", "gene_biotype") %in% colnames(annotation)
  )
  
  mr <- rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(256))
  
  stabilized_counts <- assay(vsd.obj)
  num_genes <- min(num_genes, nrow(stabilized_counts))
  stopifnot(num_genes >= 2)
  
  row_variances <- matrixStats::rowVars(stabilized_counts, useNames = TRUE)
  top_variable_genes <- stabilized_counts[order(row_variances, decreasing = TRUE)[seq_len(num_genes)], ]
  top_variable_genes <- top_variable_genes - rowMeans(top_variable_genes, na.rm = TRUE)
  
  anno_match_idx <- match(rownames(top_variable_genes), annotation$ensembl_gene_id)
  
  gene_names <- annotation$hgnc_symbol[anno_match_idx]
  missing_sym <- is.na(gene_names) | gene_names == ""
  gene_names[missing_sym] <- rownames(top_variable_genes)[missing_sym]
  gene_names <- make.unique(as.character(gene_names))
  
  biotype <- annotation$gene_biotype[anno_match_idx]
  biotype[is.na(biotype) | biotype == ""] <- "Unknown"
  row_ann <- data.frame(Biotype = factor(biotype), row.names = gene_names)
  
  rownames(top_variable_genes) <- gene_names
  
  coldata <- as.data.frame(colData(vsd.obj))
  coldata <- coldata[, !colnames(coldata) %in% c("sizeFactor", "replaceable"), drop = FALSE]
  coldata[] <- lapply(coldata, function(col) {
    if (is.numeric(col)) col else factor(as.character(col))
  })
  # Drop annotations that carry no contrast: constant, or unique per sample
  keep <- vapply(coldata, function(col) {
    n <- length(unique(col))
    n > 1 && n < ncol(vsd.obj)
  }, logical(1))
  coldata <- coldata[, keep, drop = FALSE]
  
  max_abs <- max(abs(top_variable_genes), na.rm = TRUE)
  
  pheatmap::pheatmap(
    top_variable_genes,
    color = mr,
    breaks = seq(-max_abs, max_abs, length.out = length(mr) + 1),
    annotation_col = if (ncol(coldata) > 0) coldata else NA,
    annotation_row = row_ann,
    show_colnames = TRUE,
    show_rownames = num_genes <= label_cutoff,
    fontsize_col = 6,
    fontsize_row = if (num_genes <= label_cutoff) max(4, 200 / num_genes) else 4,
    border_color = NA,
    main = title
  )
}


plot_PCA <- function(vsd.obj,
                     group_var,
                     ntop = 500,
                     label = TRUE,
                     fixed_coords = TRUE) {
  stopifnot(
    length(group_var) == 1,
    group_var %in% colnames(colData(vsd.obj))
  )
  
  ntop <- min(ntop, nrow(vsd.obj))
  pcaData <- plotPCA(vsd.obj, intgroup = group_var, ntop = ntop, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"), 1)
  
  p <- ggplot(pcaData, aes(PC1, PC2, color = .data[[group_var]])) +
    geom_point(size = 3) +
    labs(
      x = paste0("PC1: ", percentVar[1], "% variance"),
      y = paste0("PC2: ", percentVar[2], "% variance"),
      color = group_var,
      title = paste0("PCA of top ", ntop, " variable genes"),
      subtitle = paste0("Colored by ", group_var)
    )
  
  if (label) {
    p <- p + ggrepel::geom_text_repel(
      aes(label = name),
      color = "black",
      size = 2,
      max.overlaps = Inf,
      show.legend = FALSE
    )
  }
  
  if (fixed_coords) p <- p + coord_fixed()
  
  p
}