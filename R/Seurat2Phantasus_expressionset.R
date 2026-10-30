#' Aggregate expression by biological sample and cell cluster
#'
#' Sum expression within each observed sample/cluster pair while retaining
#' sparsity until the final pseudobulk matrix is constructed.
#' @param seuratobject A Seurat object.
#' @param slot Legacy alias for layer, retained for existing callers.
#' @param assay Assay name. Defaults to RNA.
#' @param sample_to_subset Biological-sample metadata column name.
#' @param cluster_to_subset Cell-cluster metadata column name.
#' @param layer Expression layer. Defaults to slot (counts by default).
#' @return A Biobase ExpressionSet with genes in rows and observed groups in
#'   columns. Annotations contain id, sample_to_subset, cluster_to_subset, n_cells.
#'   Feature annotations retain the gene identifiers in gene_id.
#' @details Use raw RNA counts for count-based downstream methods. Other layers
#'   are summed as supplied without normalization or back-transformation.
#'   Split v5 layers must first be joined explicitly with JoinLayers.
#' @export
Seurat2Phantasus_expressionset <- function(seuratobject, slot = "counts",
                                         assay = "RNA", sample_to_subset,
                                         cluster_to_subset, layer = slot) {
  check_name <- function(value, label) {
    if (!is.character(value) || length(value) != 1L || is.na(value) ||
        !nzchar(trimws(value))) {
      stop(label, " must be a single non-empty string.", call. = FALSE)
    }
  }
  if (!inherits(seuratobject, "Seurat")) {
    stop("seuratobject must be a Seurat object.", call. = FALSE)
  }
  check_name(assay, "assay")
  check_name(layer, "layer")
  check_name(sample_to_subset, "sample_to_subset")
  check_name(cluster_to_subset, "cluster_to_subset")
  if (!missing(slot) && !missing(layer) && !identical(slot, layer)) {
    stop("Specify layer or slot, or give both the same value.", call. = FALSE)
  }
  if (!assay %in% SeuratObject::Assays(seuratobject)) {
    stop("Assay not found: ", assay, call. = FALSE)
  }
  metadata <- seuratobject[[]]
  required <- c(sample_to_subset, cluster_to_subset)
  if (!all(required %in% colnames(metadata))) {
    stop("Missing metadata columns: ",
         paste(setdiff(required, colnames(metadata)), collapse = ", "),
         call. = FALSE)
  }
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    assay_object <- seuratobject[[assay]]
    layers <- SeuratObject::Layers(assay_object, search = NA)
    if (any(startsWith(layers, paste0(layer, ".")))) {
      stop("Split layers detected for '", layer,
           "'. Join the intended layers with SeuratObject::JoinLayers() first.",
           call. = FALSE)
    }
    if (!layer %in% layers) {
      stop("Layer not found: ", layer, call. = FALSE)
    }
    expression <- SeuratObject::LayerData(assay_object, layer = layer)
  } else {
    if (!layer %in% c("counts", "data", "scale.data")) {
      stop("Layer not found: ", layer,
           ". SeuratObject v4 supports counts, data, or scale.data.", call. = FALSE)
    }
    expression <- SeuratObject::GetAssayData(seuratobject, assay = assay, slot = layer)
  }
  if (length(dim(expression)) != 2L || any(dim(expression) == 0L)) {
    stop("The selected expression layer is empty.", call. = FALSE)
  }
  if (!is.matrix(expression) && !inherits(expression, "sparseMatrix")) {
    stop("The selected layer must be an in-memory matrix or sparseMatrix.",
         call. = FALSE)
  }
  cells <- colnames(expression)
  genes <- rownames(expression)
  if (is.null(cells) || anyNA(cells) || any(!nzchar(cells)) ||
      anyDuplicated(cells) || !all(cells %in% rownames(metadata))) {
    stop("Expression columns must have unique cell names matching metadata.",
         call. = FALSE)
  }
  if (is.null(genes) || anyNA(genes) || any(!nzchar(genes)) || anyDuplicated(genes)) {
    stop("Expression rows must have unique, non-empty gene names.", call. = FALSE)
  }
  # A v5 assay/layer can contain only a subset of the object's cells.
  metadata <- metadata[cells, , drop = FALSE]
  samples <- as.character(metadata[[sample_to_subset]])
  clusters <- as.character(metadata[[cluster_to_subset]])
  if (anyNA(samples) || anyNA(clusters) ||
      any(!nzchar(trimws(samples))) || any(!nzchar(trimws(clusters)))) {
    stop("Sample and cluster annotations must not contain missing or empty values.",
         call. = FALSE)
  }
  stored_values <- if (inherits(expression, "sparseMatrix")) {
    if ("x" %in% methods::slotNames(expression)) {
      methods::slot(expression, "x")
    } else {
      numeric() # Pattern sparse matrices have implicit unit entries.
    }
  } else {
    expression
  }
  if (any(!is.finite(stored_values))) {
    stop("Expression values must be finite.", call. = FALSE)
  }
  # Numeric codes keep distinct pairs separate even with punctuation in labels.
  keys <- paste(match(samples, unique(samples)),
                match(clusters, unique(clusters)), sep = ":")
  first <- which(!duplicated(keys))
  groups <- match(keys, keys[first])
  membership <- Matrix::sparseMatrix(i = seq_along(groups), j = groups,
                                   x = 1, dims = c(length(groups), length(first)))
  summed <- as.matrix(expression %*% membership)
  group_names <- make.unique(paste(clusters[first], samples[first], sep = "."))
  dimnames(summed) <- list(genes, group_names)
  annotations <- data.frame(id = group_names, sample_to_subset = samples[first],
                            cluster_to_subset = clusters[first],
                            n_cells = tabulate(groups, nbins = length(first)),
                            row.names = group_names, stringsAsFactors = FALSE)
  Biobase::ExpressionSet(assayData = summed,
                        phenoData = methods::new("AnnotatedDataFrame", annotations),
                        featureData = methods::new("AnnotatedDataFrame",
                          data.frame(gene_id = genes, row.names = genes,
                                     stringsAsFactors = FALSE)))
}
