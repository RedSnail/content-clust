test_i = c(2, 3, 4, 5, 6, 5, 6, 2, 3, 3, 4, 6, 0, 1)
test_p = c(0, 5, 7, 9, 11, 12, 14)
gene_names =  c("a", "a.a", "shit1", "shit2", "a.a.a", "b")
expression = test_p[-1] - test_p[-length(test_p)]
raw = ContentClust(test_i, test_p, 7, gene_names, order(test_p[-1] - test_p[-length(test_p)], decreasing = T) - 1)
raw

gene_set = readRDS("~/R/content_clust/data/vargenes_indices.rds")
var_genes_mat = data[gene_set,]
to_content_clust = t(var_genes_mat)
cells_expr = to_content_clust@p[-1] - to_content_clust@p[-length(to_content_clust@p)]
cells_clust_raw = ContentClust(to_content_clust@i, to_content_clust@p, dim(to_content_clust)[1],
                               rownames(to_content_clust), 
                               order(cells_expr, decreasing = T) - 1)


raw_to_json = function(cdata) {
  out = "{%s}"
  
  
}
