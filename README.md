# content-clust

here we try different approaches in scRNAseq that have one idea in common: zero gene expression has qualitive difference from any other value.\

### current directions of work:
######1. gene hierarchy
Some genes are expressed only in case another gene is expressed. So we make an orgraph of genes where gene A is connected to gene B if all cells that express B, express A. Then we analyse the graph. You can see the full algorithm in src/source.cpp
######2. Gower distance
Gower distance is a measure of simmilarity between objects by both quantitive and qualitive features. By now, no code on this theme is uploaded.



