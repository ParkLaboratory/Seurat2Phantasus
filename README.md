# About Seurat2Phantasus
The tool performs a pseudobulk conversion of the single cell data of a Seurat object into an expression set object or gct file suitable for the Phantasus tool (https://artyomovlab.wustl.edu/phantasus/).
All this allows to perform GUI-based pseudobulk analysis of single-cell RNA sequencing data using all the advantages of the Phantasus tool, including PCA analysis, Differential expression (DESeq, Limma), 
Pathways analysis (fgsea), and others (for the complete list of functions and manual, see Phantasus help: https://ctlab.github.io/phantasus-doc/quick-start.html). 
Another helpful feature of this package is the ability to generate an expression set object compatible with many bulk RNA seq analysis tools, the example below will show how to remove the batch effect using Combat-Seq.

### Install from GitHub
#### Install dependencies:
```r
install.packages('Seurat')
BiocManager::install(c('a4Base','Biobase','phantasus'))
```
#### Install Seurat2Phantasus
```r
install.packages("devtools")
devtools::install_github("ParkLaboratory/Seurat2Phantasus")
```
## Usage

Basic usage (users need to input two cells annotation (for example: sample and cluster), indicate slot and assay, and output file location):

```r
Seurat2Phantasus_gct(seuratobject = MySeuratObject,
           slot = "counts",
           assay = "RNA",
           sample_to_subset = "sample",
           cluster_to_subset = "seurat_clusters",
           output_gct_file = "MyPseudobulk.gct")
```
The resulting .gct file could be uploaded to Phantasus tool (https://artyomovlab.wustl.edu/phantasus/) for future analysis.

Expression set generation (for advanced analysis).

```r
es <- Seurat2Phantasus_expressionset(seuratobject = MySeuratObject,
            slot = "counts",
            assay = "RNA",
            sample_to_subset = "sample",
            cluster_to_subset = "seurat_clusters")

#add information about batch
pData(es)["batch"] <- c("batch1","batch1","batch1","batch2","batch2","batch2")

#remove batch effect using ComBat_seq (https://github.com/zhangyuqing/ComBat-seq)
export <- ComBat_seq(exprs(es), batch = pData(es)[,"batch", drop=TRUE], group = NULL)

#generate batch-free .gct
phantasus::write.gct(Biobase::ExpressionSet(assayData = export, phenoData = new("AnnotatedDataFrame",pData(es))), 
          "MyPseudobulk_wo_batcheffect.gct", gzip = FALSE)

```
