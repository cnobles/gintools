context("Generate position IDs")

chr <- c("chr1", "chr3", "chrX")
strands <- c("+", "-", "+")
starts <- c(900231, 13254892, 603292)
ends <- c(900431, 13255292, 603592)
ranges <- IRanges::IRanges(start = starts, end = ends)
input_gr <- GenomicRanges::GRanges(
  seqnames = chr, ranges = ranges, strand = strands
)

expected_posids <- c("chr1+900231", "chr3-13255292", "chrX+603292")

test_that(
  desc = "Output of posids from a GRange",
  code = {

    output <- generate_posid(input_gr)
    expect_equal(output, expected_posids)

  }
)
