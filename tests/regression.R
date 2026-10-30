library(Seurat2Phantasus)

expect_error <- function(code, text) {
  message <- tryCatch({ force(code); NULL }, error = conditionMessage)
  stopifnot(!is.null(message), grepl(text, message, fixed = TRUE))
}
counts <- matrix(c(1, 0, 2, 3, 1, 0, 0, 4, 2, 2, 1, 1, 5, 2, 0, 1, 1, 3),
                 nrow = 3, dimnames = list(paste0("g", 1:3), paste0("c", 1:6)))
metadata <- data.frame(sample = c("s1", "s1", "s2", "s2", "s2", "s1"),
                       cluster = c("T", "T", "T", "B", "T", "B"),
                       row.names = colnames(counts))
make_object <- function(x = counts, meta = metadata) {
  SeuratObject::CreateSeuratObject(counts = x, meta.data = meta)
}
aggregate_es <- function(object, ...) {
  Seurat2Phantasus_expressionset(object, sample_to_subset = "sample",
                                cluster_to_subset = "cluster", ...)
}
expected <- cbind(c(4, 1, 2), c(5, 6, 2), c(2, 1, 1), c(1, 1, 3))
dimnames(expected) <- list(rownames(counts), c("T.s1", "T.s2", "B.s2", "B.s1"))
object <- make_object()
es <- aggregate_es(object)
stopifnot(identical(Biobase::exprs(es), expected),
          identical(Biobase::fData(es)$gene_id, rownames(counts)),
          identical(Biobase::pData(es)$id, colnames(expected)),
          identical(as.integer(Biobase::pData(es)$n_cells), c(2L, 2L, 1L, 1L)),
          identical(Biobase::pData(es)$sample_to_subset, c("s1", "s2", "s2", "s1")),
          identical(Biobase::exprs(aggregate_es(object, slot = "counts")), expected))
sparse_object <- make_object(Matrix::Matrix(counts, sparse = TRUE))
stopifnot(identical(Biobase::exprs(aggregate_es(sparse_object)), expected))
one <- make_object(counts[, 1, drop = FALSE], metadata[1, , drop = FALSE])
stopifnot(identical(unname(Biobase::exprs(aggregate_es(one))),
                    unname(counts[, 1, drop = FALSE])))

# Labels that produce the same printable name must still form distinct groups.
collision <- metadata
collision$sample <- c("b.c", "c", "x", "x", "x", "x")
collision$cluster <- c("a", "a.b", "z", "z", "z", "z")
collision_es <- aggregate_es(make_object(meta = collision))
stopifnot(ncol(collision_es) == 3L,
          !anyDuplicated(colnames(collision_es)),
          identical(unname(Biobase::exprs(collision_es)[, 1:2]),
                    unname(counts[, 1:2])))

expect_error(aggregate_es(list()), "must be a Seurat object")
expect_error(aggregate_es(object, assay = "missing"), "Assay not found")
expect_error(aggregate_es(object, layer = "missing"), "Layer not found")
expect_error(aggregate_es(object, slot = "data", layer = "counts"), "Specify layer or slot")
bad <- object
bad[["sample"]] <- NULL
expect_error(aggregate_es(bad), "Missing metadata columns")
bad <- object
bad$sample <- c(NA, metadata$sample[-1])
expect_error(aggregate_es(bad), "missing or empty")
bad$sample <- c(" ", metadata$sample[-1])
expect_error(aggregate_es(bad), "missing or empty")
expect_error(Seurat2Phantasus_gct(object, sample_to_subset = "sample",
                                cluster_to_subset = "cluster", output_gct_file = ""),
             "non-empty path")

output <- tempfile(fileext = ".gct")
result <- Seurat2Phantasus_gct(object, sample_to_subset = "sample",
                             cluster_to_subset = "cluster", output_gct_file = output)
stopifnot(identical(result, output), file.exists(output),
          identical(readLines(output, n = 1L), "#1.3"))
roundtrip <- phantasus::read.gct(output)
stopifnot(isTRUE(all.equal(Biobase::exprs(roundtrip), expected,
                           check.attributes = TRUE)),
          identical(as.character(Biobase::pData(roundtrip)$sample_to_subset),
                    Biobase::pData(es)$sample_to_subset))
unlink(output)

if (utils::packageVersion("SeuratObject") >= "5.0.0") {
  split_assay <- SeuratObject::CreateAssay5Object(
    counts = list(a = counts[, 1:3], b = counts[, 4:6]))
  split_object <- SeuratObject::CreateSeuratObject(counts = split_assay,
                                                meta.data = metadata)
  expect_error(aggregate_es(split_object), "Split layers detected")
  joined <- SeuratObject::JoinLayers(split_object, assay = "RNA")
  stopifnot(identical(Biobase::exprs(aggregate_es(joined)), expected))
  subset_object <- object
  subset_object[["subset"]] <- SeuratObject::CreateAssay5Object(
    counts = counts[, c(1, 3, 5), drop = FALSE])
  subset_es <- aggregate_es(subset_object, assay = "subset")
  stopifnot(ncol(subset_es) == 2L,
            identical(as.integer(Biobase::pData(subset_es)$n_cells), c(1L, 2L)),
            identical(unname(Biobase::exprs(subset_es)),
                      unname(cbind(c(1, 0, 2), c(5, 6, 2)))))
}
# Run the actual README example in an isolated temporary directory.
if (file.exists("README.md")) {
  readme <- readLines("README.md")
  start <- which(readme == "```r")[2L]
  end <- which(readme == "```" & seq_along(readme) > start)[1L]
  stopifnot(!is.na(start), !is.na(end))
  example_code <- readme[seq.int(start + 1L, end - 1L)]
  local({
    example_dir <- tempfile("seurat2phantasus-example-")
    dir.create(example_dir)
    original_dir <- setwd(example_dir)
    on.exit({ setwd(original_dir); unlink(example_dir, recursive = TRUE) })
    example_env <- new.env()
    eval(parse(text = example_code), envir = example_env)
    stopifnot(file.exists("pseudobulk.gct"), ncol(example_env$es) == 3L)
  })
}
message("Pseudobulk regression checks passed.")
