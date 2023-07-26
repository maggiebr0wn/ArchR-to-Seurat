<img src = "https://github.com/maggiebr0wn/ArchR-to-Seurat/blob/main/archR2seurat.png">


A tool which converts an ArchR object to Seurat object.

To run:

    Rscript ./archR2seurat.R <archr_project_path>

Written for an ArchR object meeting the following criteria:

<li> single cell Multiomic data: scRNA-seq and scATAC-seq from the same cells. </li>
<li> The two matrices to be transferred are "PeakMatrix" and "GeneExpressionMatrix".</li>
<li> The UMAP embedding from ArchR is called "UMAP". </li>
<p>
<br>
These criteria may not apply to all ArchR objects; matrices of other types and names can be transferred as well using the getMatrixFromProject() function. Same with UMAP embeddings of other names, using the getEmbedding() function.
</p>
<br>
Cheers!
