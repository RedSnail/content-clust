test_i <- c(2, 3, 4, 5, 6, 5, 6, 2, 3, 3, 4, 6, 0, 1)
#           0  1  2  3  4  5  6  7  8  9 10 11 12 13
test_p <- c(0, 5, 7, 9, 11, 12, 14)
#           0  1  2  3  4   5   6
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
	expect_equal(substract_sorted(c(1, 4, 5), test_i, c(12), c(14)), c(4, 5))
	expect_equal(substract_sorted(c(2, 3, 4, 5, 6), test_i, c(5), c(7)), c(2, 3, 4))
	expect_equal(substract_sorted(c(2, 3, 4, 5, 6), test_i, c(5, 2), c(7, 3)), c(2, 3)) 
})

test_that("test other simple vector functions", {
  # test is_in()
  expect_equal(is_in(1, c(1, 2, 3)), TRUE)
  expect_equal(is_in(4, c(1, 2, 3)), FALSE)
  
  # test concat()
  expect_equal(concat(c("a", "b", "c"), " "), "a b c")
  expect_equal(concat(c("aa", "bb", "cc"), ""), "aabbcc")
})

test_that("count each gene expression", {
	# test expression()
  expect_equal(expression(test_p), c(5, 2, 2, 2, 1, 2))

})

gene_order <- order(expression(test_p), decreasing=TRUE) - 1

test_that("suborgraph processing preparations", {
	# test children finding
	expect_equal(find_children(c(-1, 0, 1, -1, 3), c(-1)), c(0, 3))
	expect_equal(find_children(c(-1, 0, 1, -1, 3), c(0)), c(1))
	expect_equal(find_children(c(-1, 0, 1, -1, 3), c(2)), vector(mode="numeric", length=0))
	expect_equal(find_children(c(-1, 0, 1, -1, 3), c(-1, 0)), c(0, 1, 3))
	expect_equal(find_children(c(-1, 0, 1, -1, 3), c(0, 2)), c(1))


	# test cell set initiation
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, c(-1)), c(0, 1, 2, 3, 4, 5, 6))
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, c(5)), c(6))
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, c(0)), c(2, 3, 4, 5, 6))
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, c(-1, 0)), c(0, 1, 2, 3, 4, 5, 6))
	expect_equal(get_cell_set(test_i, test_p, 7, gene_order, c(2, 3)), c(2, 3, 4))
})

test_that("genes that have no overlaps are selected correctly", {
	# test getting sep genes
	expect_equal(get_sep_genes(test_i, test_p, gene_order, c(5, 2, 3)), list(c(5), c(2, 3)))
})

test_that("whole algoritm works okay", {
	gene_names <- c("a", "a.a", "shit1", "shit2", "a.a.a", "b")
	# test subgraph processing
	outlist <- process_oriented_graph(test_i, test_p, 7, gene_names, gene_order, c(-1, 0, 0, 0, -1, 1), c(0))
	expected_outlist <-list("a.a" = list("a.a.a" = c(6), "unclassified" = c(5)), "shit1 shit2" = c(2, 3, 4))
	expect_equal(outlist, expected_outlist)
	whole_p <- c(0, 3, 5, 7, 9, 10, 12)
  whole_i <- c(2, 3, 4, 4, 5, 6, 7, 7, 8, 7, 0, 1)
  whole_order <- order(expression(whole_p), decreasing = TRUE) - 1

  whole_names <- c("a.1", "a.2", "b.1", "b.2", "b.a.1", "c.1")


	# finally, test whole algoritm
	outlist_2 <- content_clust(whole_i, whole_p, 9, whole_names, whole_order)
	expected_outlist_2 <- list("a.1 a.2" = c(2, 3, 4, 5), "b.1 b.2" = list("b.a.1"=c(7), "unclassified"=c(6, 8)), "c.1" = c(0, 1))
	expect_equal(outlist_2, expected_outlist_2)
})
