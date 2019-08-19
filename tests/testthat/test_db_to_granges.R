context("Database structure to GRanges object")

input_df <- data.frame(
  "chr" = c("chr1", "chr2", "chr2", "chr3"),
  "position" = c(5379927, 92775920, 2719573, 7195924),
  "breakpoint" = c(5380070, 92775995, 2719450, 7195890),
  "strand" = c("+", "+", "-", "-"),
  "sampleName" = c("GTSP0001-3", "GTSP0003-1", "GTSP0002-4", "GTSP0003-1"),
  "is.multihit" = c(FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)

output_gr <- GenomicRanges::GRanges(
  seqnames = factor(c("chr1", "chr2", "chr2", "chr3")),
  ranges = IRanges::IRanges(
    start = c(5379927, 92775920, 2719450, 7195890),
    end = c(5380070, 92775995, 2719573, 7195924)
  ),
  strand = factor(c("+", "+", "-", "-"), levels = c("+", "-", "*")),
  samplename = c("GTSP0001-3", "GTSP0003-1", "GTSP0002-4", "GTSP0003-1"),
  specimen = c("GTSP0001", "GTSP0003", "GTSP0002", "GTSP0003"),
  is.multihit = c(FALSE, TRUE, FALSE, TRUE)
)

test_that("Converted database dataframe to GRanges object", {
  
  set.seed(1)
  res <- db_to_granges(input_df)
  expect_equal(
    as.data.frame(res), 
    as.data.frame(output_gr)
  )
  
})
