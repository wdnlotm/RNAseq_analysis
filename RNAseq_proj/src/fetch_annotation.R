# Fetch Ensembl gene annotation (ID -> symbol -> biotype) and cache it to CSV.
#
# Why this does not use biomaRt:
#   Ensembl's new website retired the classic BioMart path on www.ensembl.org
#   (it returns {"status_code": 404, "details": "No supported new Ensembl
#   equivalent for this URL"}), and the regional mirrors useast/asia now return
#   403. biomaRt 2.62's .getArchiveList() hardcodes exactly those three hosts, so
#   useEnsembl()/useMart() fail with "Unable to contact any Ensembl mirror"
#   regardless of the `host` argument. The BioMart *service* itself is alive on
#   the dated archive hosts, so we speak its TSV/XML protocol directly.
#
# Pinning to a dated archive is also better for reproducibility than tracking
# whatever "current" resolves to.

ENSEMBL_ARCHIVE <- "https://may2025.archive.ensembl.org"  # Ensembl 114, GRCh38
BIOMART_DATASET <- "hsapiens_gene_ensembl"                # mouse: mmusculus_gene_ensembl
BIOMART_SYMBOL <- "hgnc_symbol"                           # mouse: mgi_symbol

#' Download a gene annotation table from an Ensembl BioMart archive
#'
#' Retrieves every gene in the dataset rather than filtering on a list of IDs,
#' which keeps the request to a single call and lets callers join locally.
#'
#' @param archive Base URL of a dated Ensembl archive host.
#' @param dataset BioMart dataset name.
#' @param symbol_attribute Species-appropriate symbol attribute.
#' @param timeout_sec Request timeout; the full human table is a few MB.
#'
#' @return A data frame with `ensembl_gene_id`, `hgnc_symbol`, `gene_biotype`,
#'   one row per gene. Genes with no symbol get their Ensembl ID as the symbol,
#'   matching what the downstream plotting and enrichment code expects.
fetch_ensembl_annotation <- function(archive = ENSEMBL_ARCHIVE,
                                    dataset = BIOMART_DATASET,
                                    symbol_attribute = BIOMART_SYMBOL,
                                    timeout_sec = 300) {
  query_xml <- paste0(
    '<?xml version="1.0" encoding="UTF-8"?><!DOCTYPE Query>',
    '<Query virtualSchemaName="default" formatter="TSV" header="1" ',
    'uniqueRows="1" count="" datasetConfigVersion="0.6">',
    '<Dataset name="', dataset, '" interface="default">',
    '<Attribute name="ensembl_gene_id"/>',
    '<Attribute name="', symbol_attribute, '"/>',
    '<Attribute name="gene_biotype"/>',
    '</Dataset></Query>'
  )

  resp <- httr::POST(
    paste0(archive, "/biomart/martservice"),
    body = list(query = query_xml),
    encode = "form",
    httr::timeout(timeout_sec)
  )

  if (httr::status_code(resp) != 200L) {
    stop("BioMart request to ", archive, " failed with status ",
         httr::status_code(resp), call. = FALSE)
  }

  txt <- httr::content(resp, as = "text", encoding = "UTF-8")

  # BioMart reports query errors in a 200 response body
  if (grepl("Query ERROR", txt, fixed = TRUE)) {
    stop("BioMart returned an error: ", substr(txt, 1, 300), call. = FALSE)
  }

  out <- readr::read_tsv(I(txt), show_col_types = TRUE, progress = FALSE)
  if (ncol(out) != 3L) {
    stop("Expected 3 columns from BioMart, got ", ncol(out), call. = FALSE)
  }
  names(out) <- c("ensembl_gene_id", "hgnc_symbol", "gene_biotype")

  out <- dplyr::distinct(out, ensembl_gene_id, .keep_all = TRUE)

  missing_symbol <- is.na(out$hgnc_symbol) | out$hgnc_symbol == ""
  out$hgnc_symbol[missing_symbol] <- out$ensembl_gene_id[missing_symbol]

  stopifnot(!any(duplicated(out$ensembl_gene_id)))
  out
}

#' Load the cached annotation, downloading it on first use
#'
#' @param path CSV cache location.
#' @param refresh Force a re-download even if the cache exists.
load_annotation <- function(path, refresh = FALSE, ...) {
  if (!refresh && file.exists(path)) {
    return(read.csv(path, header = TRUE, stringsAsFactors = FALSE))
  }
  annotation <- fetch_ensembl_annotation(...)
  dir.create(dirname(path), showWarnings = FALSE, recursive = TRUE)
  readr::write_csv(annotation, path)
  message("Wrote ", nrow(annotation), " genes to ", path)
  annotation
}
