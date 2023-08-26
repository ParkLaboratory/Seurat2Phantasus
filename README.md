# About Seurat2Phantasus
The tool performs a pseudobulk conversion of the single cell data of a Seurat object into an expression set object or gct file suitable for the Phantasus tool (https://artyomovlab.wustl.edu/phantasus/).
All this allows to perform GUI-based pseudobulk analysis of single-cell RNA sequencing data using all the advantages of the Phantasus tool, including PCA analysis, Differential expression (DESeq, Limma), 
Pathways analysis (fgsea), and others (for the complete list of functions and manual see Phantasus help: https://ctlab.github.io/phantasus-doc/quick-start.html). Another helpful feature of this package is the ability to generate an expression set object compatible with many bulk RNA seq analysis tools, the example below will show how to remove the donor effect using Combat-Seq (https://github.com/zhangyuqing/ComBat-seq).

### Install from GitHub
#### Install dependencies
devtools::install_github("ParkLaboratory/Seurat2Phantasus")
