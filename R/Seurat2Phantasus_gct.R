#' Export sample/cluster pseudobulk profiles to a GCT file
#' @inheritParams Seurat2Phantasus_expressionset
#' @param output_gct_file Output path. Its directory must exist.
#' @return Invisibly returns the output path after writing GCT version 1.3.
#' @export
Seurat2Phantasus_gct <- function(seuratobject, slot = "counts", assay = "RNA",
                                sample_to_subset, cluster_to_subset,
                                output_gct_file, layer = slot) {
  if (!is.character(output_gct_file) || length(output_gct_file) != 1L ||
      is.na(output_gct_file) || !nzchar(trimws(output_gct_file))) {
    stop("output_gct_file must be a single non-empty path.", call. = FALSE)
  }
  if (!dir.exists(dirname(output_gct_file))) {
    stop("Output directory does not exist: ", dirname(output_gct_file), call. = FALSE)
  }
  if (!missing(slot) && !missing(layer) && !identical(slot, layer)) {
    stop("Specify layer or slot, or give both the same value.", call. = FALSE)
  }
  es <- Seurat2Phantasus_expressionset(
    seuratobject = seuratobject, assay = assay,
    sample_to_subset = sample_to_subset, cluster_to_subset = cluster_to_subset,
    layer = layer)
  phantasusLite::writeGct(es, output_gct_file, gzip = FALSE)
  invisible(output_gct_file)
}
