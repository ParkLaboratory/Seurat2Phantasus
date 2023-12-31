% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Seurat2Phantasus_gct.R
\name{Seurat2Phantasus_gct}
\alias{Seurat2Phantasus_gct}
\title{A Seurat2Phantasus_gct Function
This function allows you to make pseudobulk gct file for Phantasus}
\usage{
Seurat2Phantasus_gct(
  seuratobject,
  slot = "counts",
  assay = "RNA",
  sample_to_subset,
  cluster_to_subset,
  output_gct_file
)
}
\arguments{
\item{seuratobject}{Your seuratobject}

\item{slot}{Slot: counts - raw, data - log2 Default is counts slot}

\item{assay}{Assay: RNA, SCT, etc Default is RNA slot}

\item{sample_to_subset}{sample metadata annotation}

\item{cluster_to_subset}{cluster metadata annotation}

\item{output_gct_file}{output gct file}
}
\description{
A Seurat2Phantasus_gct Function
This function allows you to make pseudobulk gct file for Phantasus
}
