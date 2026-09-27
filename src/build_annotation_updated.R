# Build an annotation table keyed to the GTF used for quantification.
#
# The RSEM merged counts file carries the gene_id and gene_name straight from the
# GTF that FASTQs were processed against, so it covers every gene in the count
# matrix (58,884/58,884) including ~1,290 IDs that have since been retired from
# Ensembl. It has no biotype column, so biotypes are joined in from the cached
# Ensembl BioMart table. Retired IDs therefore get a symbol but a missing
# biotype, which is the honest representation: they have no current biotype.
#
# Symbols are NOT unique (~1,500 duplicated gene_name values in the GTF), so
# ensembl_gene_id remains the key for all joins and rownames.

#' Build the GTF-anchored annotation table
#'
#' @param rsem_file RSEM merged gene counts TSV; only the first two columns
#'   (gene_id, gene_name) are read.
#' @param biotype_file Cached Ensembl annotation CSV supplying `gene_biotype`.
#'
#' @return A data frame with one row per gene in the GTF:
#'   `ensembl_gene_id` (version stripped, the join key), `gene_id` (as written in
#'   the GTF), `hgnc_symbol`, `gene_biotype` (NA for retired IDs).
build_annotation_updated <- function(rsem_file, biotype_file) {
  gtf <- readr::read_tsv(rsem_file, col_select = 1:2, show_col_types = FALSE)
  stopifnot(identical(names(gtf), c("gene_id", "gene_name")))

  gtf$ensembl_gene_id <- sub("\\..*$", "", gtf$gene_id)

  # Version stripping collapses a handful of IDs; the duplicate rows are
  # identical in the counts, so keeping the first is safe.
  gtf <- dplyr::distinct(gtf, ensembl_gene_id, .keep_all = TRUE)

  biotypes <- read.csv(biotype_file, stringsAsFactors = FALSE)
  stopifnot(all(c("ensembl_gene_id", "gene_biotype") %in% names(biotypes)))
  biotypes <- dplyr::distinct(biotypes, ensembl_gene_id, .keep_all = TRUE)

  out <- gtf |>
    dplyr::left_join(
      dplyr::select(biotypes, ensembl_gene_id, gene_biotype),
      by = "ensembl_gene_id"
    ) |>
    dplyr::mutate(
      # A GTF gene_name is occasionally blank; fall back to the ID so plotting
      # and enrichment code always has a usable label.
      hgnc_symbol = dplyr::if_else(
        is.na(gene_name) | gene_name == "", ensembl_gene_id, gene_name
      )
    ) |>
    dplyr::select(ensembl_gene_id, gene_id, hgnc_symbol, gene_biotype)

  stopifnot(!any(duplicated(out$ensembl_gene_id)))
  out
}

#' Write the table to CSV, rebuilding only when missing or when forced
load_annotation_updated <- function(path, rsem_file, biotype_file,
                                    refresh = FALSE) {
  if (!refresh && file.exists(path)) {
    return(read.csv(path, stringsAsFactors = FALSE))
  }
  out <- build_annotation_updated(rsem_file, biotype_file)
  readr::write_csv(out, path)
  message("Wrote ", nrow(out), " genes to ", path, " (",
          sum(is.na(out$gene_biotype)), " without a biotype)")
  out
}
