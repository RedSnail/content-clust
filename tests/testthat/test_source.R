test_i <- c(2, 3, 4, 5, 6, 5, 6, 2, 3, 3, 4, 6, 0, 1)
test_p <- c(0, 5, 7, 9, 11, 12, 14)
# gene_names =  c("a", "a.a", "shit1", "shit2", "a.a.a", "b")
# expression = test_p[-1] - test_p[-length(test_p)]
# raw = ContentClust(test_i, test_p, 7, gene_names, order(test_p[-1] - test_p[-length(test_p)], decreasing = T) - 1)

test_that("functions that work with sets like sorted vectors", {
	# test small_in_big()
	expect_equal(small_in_big(test_i, 5, 7, 0, 5), TRUE)
	expect_equal(small_in_big(test_i, 0, 5, 5, 7), FALSE)
	expect_equal(small_in_big(test_i, 12, 14, 0, 5), FALSE)

	# test have_overlap()
	expect_equal(have_overlap(test_i, 5, 7, 7, 9), FALSE)
	expect_equal(have_overlap(test_i, 7, 9, 9, 11), TRUE)

	# test substract_sorted()
	expect_equal(substract_sorted(c(1, 4, 5), test_i, 12, 14), c(4, 5)) 
})

test_that("count each gene expression", {
	# test expression()
        expect_equal(expression(test_p), c(5, 2, 2, 2, 1, 2))

})

gene_order <- order(expression(test_p), decreasing=TRUE) - 1

test_that("suborgraph processing preparations", {
	# test children finding
	expect_equal(find_children(c(-1, 0, 1, -1, 3), -1), c(0, 3))
	expect_equal(find_children(c(-1, 0, 1, -1, 3), 0), c(1))
	expect_equal(find_children(c(-1, 0, 1, -1, 3), 2), vector(mode="numeric", length=0))


	# test cell set initiation
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, -1), c(0, 1, 2, 3, 4, 5, 6))
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, 5), c(6))
})

test_that("genes that have no overlaps are selected correctly", {
	# test getting sep genes
	expect_equal(get_sep_genes(test_i, test_p, gene_order, c(5, 2, 3)), c(5))
})

test_that("whole algoritm works okay", {
	gene_names <- c("a", "a.a", "shit1", "shit2", "a.a.a", "b")
	# test subgraph processing
	outlist <- process_oriented_graph(test_i, test_p, 7,
					 gene_names, gene_order, 
					 c(-1, 0, 0, 0, -1, 1), 0)
	expected_outlist <-list("a.a" = list("a.a.a" = c(6), "unclussified" = c(5)), "unclussified" = c(2, 3, 4))
	expect_equal(outlist, expected_outlist)

	# finally, test whole algoritm
	outlist_2 <- content_clust(test_i, test_p, 7, gene_names, gene_order)
	expected_outlist_2 <- list("a"=outlist, "b"=c(0, 1))
	expect_equal(outlist_2, expected_outlist_2)
})
