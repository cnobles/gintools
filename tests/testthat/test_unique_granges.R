context("Unique GRanges")

input_gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(factor("chr1"), 20),
  ranges = IRanges::IRanges(
    start = c(
      598L, 598L, 599L, 599L, 599L, 599L, 599L, 599L, 600L, 600L, 
      600L, 600L, 600L, 600L, 600L, 600L, 600L, 600L, 600L, 600L
    ),
    end = c(
      687L, 687L, 687L, 687L, 687L, 687L, 687L, 687L, 660L, 660L, 
      640L, 640L, 640L, 640L, 640L, 640L, 640L, 640L, 640L, 640L
    )
  ),
  strand = S4Vectors::Rle(
    factor(c("+", "-"), levels = c("+", "-", "*")), 
    lengths = 10
  ),
  sample = rep(LETTERS[1:2], 10),
  counts = 1
)

output_gr <- GenomicRanges::GRanges(
  seqnames = factor("chr1"),
  ranges = IRanges::IRanges(
    start = c(598L, 598L, 599L, 599L, 600L, 600L, 600L, 600L),
    end = c(687L, 687L, 687L, 687L, 640L, 640L, 660L, 660L)
  ),
  strand = factor(
    c("+", "+", "+", "+", "-", "-", "+", "+"), 
    levels = c("+", "-", "*")
  ),
  sample = c("A", "B", "A", "B", "A", "B", "A", "B"),
  counts = c(1, 1, 3, 3, 5, 5, 1, 1)
)

test_that("Unique GRanges expected in output", {
  
  set.seed(1)
  res <- unique_granges(input_gr, sum.cols = "counts")
  
  expect_equal(
    as.data.frame(res, row.names = NULL),
    as.data.frame(output_gr, row.names = NULL)
  )

})
