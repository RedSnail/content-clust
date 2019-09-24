test_i = c(2, 3, 4, 5, 6, 5, 6, 2, 3, 3, 4, 6, 0, 1)
test_p = c(0, 5, 7, 9, 11, 12, 14)
# gene_names =  c("a", "a.a", "shit1", "shit2", "a.a.a", "b")
# expression = test_p[-1] - test_p[-length(test_p)]
# raw = ContentClust(test_i, test_p, 7, gene_names, order(test_p[-1] - test_p[-length(test_p)], decreasing = T) - 1)

test_that("auxilary functions work okay", {
	# test expression()
	expect_equal(expression(test_p), c(5, 2, 2, 2, 1, 2))

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

