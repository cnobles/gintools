context("Refine breakpoints")

input_gr <- GenomicRanges::GRanges(
  seqnames = S4Vectors::Rle(factor("chr1"), 40),
  ranges = IRanges::IRanges(
    start = c(
      598L, 598L, 599L, 599L, 599L, 599L, 599L, 599L, 600L, 600L, 
      600L, 600L, 600L, 600L, 600L, 600L, 601L, 601L, 601L, 602L, 
      603L, 603L, 603L, 604L, 604L, 615L, 615L, 615L, 616L, 617L, 
      575L, 575L, 577L, 577L, 579L, 573L, 573L, 573L, 572L, 574L
    ),
    end = c(
      686L, 687L, 687L, 687L, 686L, 658L, 658L, 658L, 658L, 657L, 
      639L, 639L, 639L, 638L, 637L, 676L, 677L, 676L, 675L, 675L, 
      657L, 657L, 658L, 658L, 658L, 659L, 659L, 659L, 659L, 659L, 
      659L, 659L, 660L, 660L, 660L, 660L, 661L, 661L, 661L, 661L
    )
  ),
  strand = S4Vectors::Rle(
    factor(c("+", "-"), levels = c("+", "-", "*")), 
    lengths = 20
  )
)


test_that("Standardized end positions are expected", {
  
  set.seed(1)
  res <- refine_breakpoints(input_gr)
  
  calculated_bps <- GenomicRanges::start(
    GenomicRanges::flank(res, -1, start = FALSE)
  )
  
  expected_bps <- rep(
    c(687L, 658L, 639L, 676L, 603L, 615L, 575L, 573L), each = 5
  )
  
  expect_equal(calculated_bps, expected_bps)
  
})
