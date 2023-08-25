#' A Seurat2Phantasus_expressionset Function
#' This function allows you to make pseudobulk expression set file
#' @param seuratobject Your seuratobject
#' @param slot Slot: counts - raw, data - log2 Default is counts slot
#' @param assay Assay: RNA, SCT, etc Default is RNA slot
#' @param sample_to_subset sample metadata annotation
#' @param cluster_to_subset cluster metadata annotation
#' @keywords
#' @examples
#' @export

Seurat2Phantasus_expressionset <-
  function(seuratobject,
           slot = "counts",
           assay = "RNA",
           sample_to_subset,
           cluster_to_subset)

  {
    SeuratObject.es = Biobase::ExpressionSet(assayData = as.matrix(Seurat::GetAssayData(seuratobject, slot = slot, assay = assay)), phenoData =  new("AnnotatedDataFrame",seuratobject@meta.data))

    sample_to_subset.all.es <- Biobase::ExpressionSet(assayData = as.matrix(Biobase::fData(SeuratObject.es)))

    for (sample in unique(SeuratObject.es[[sample_to_subset]])) {

      sample_to_subset.es <- SeuratObject.es[,SeuratObject.es[[sample_to_subset]] == sample]

      cluster_to_subset.all.es <- Biobase::ExpressionSet(assayData = as.matrix(Biobase::fData(SeuratObject.es)))

      for (cluster in unique(sample_to_subset.es[[cluster_to_subset]])) {

        cluster_to_subset.es <- sample_to_subset.es[ , sample_to_subset.es[[cluster_to_subset]] == cluster]

        print(paste0("Sample: ",sample," cluster: ",cluster," contains: ",ncol(exprs(cluster_to_subset.es))," cells"))

        cluster_to_subset.es.rowsum <- Biobase::ExpressionSet(assayData = as.matrix(rowSums(exprs(cluster_to_subset.es))))

        colnames(cluster_to_subset.es.rowsum) <- paste0(cluster,'.',sample)
        pData(cluster_to_subset.es.rowsum)$sample_to_subset = sample
        pData(cluster_to_subset.es.rowsum)$cluster_to_subset = cluster

        cluster_to_subset.all.es <- a4Base::combineTwoExpressionSet(cluster_to_subset.es.rowsum,cluster_to_subset.all.es)
        rm(cluster_to_subset.es.rowsum)
      }
      sample_to_subset.all.es <- a4Base::combineTwoExpressionSet(cluster_to_subset.all.es,sample_to_subset.all.es)
      rm(cluster_to_subset.all.es)
    }
    sample_to_subset.all.es
  }
