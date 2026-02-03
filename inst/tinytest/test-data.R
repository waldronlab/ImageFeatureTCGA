# Test TCGAcodesAvailable dataset

# Load data
data("TCGAcodesAvailable", package = "ImageFeatureTCGA")

# Test class
expect_inherits(TCGAcodesAvailable, "data.frame")

# Test dimensions
expect_equal(nrow(TCGAcodesAvailable), 31)
expect_equal(ncol(TCGAcodesAvailable), 4)

# Test column names
expected_cols <- c("diseaseCodes", "slide_level_available", 
                   "tile_level_available", "hover_available")
expect_true(all(expected_cols %in% names(TCGAcodesAvailable)))

# Test column types
expect_inherits(TCGAcodesAvailable$diseaseCodes, "character")
expect_inherits(TCGAcodesAvailable$slide_level_available, "logical")
expect_inherits(TCGAcodesAvailable$tile_level_available, "logical")
expect_inherits(TCGAcodesAvailable$hover_available, "logical")

# Test disease codes format (all start with TCGA_)
expect_true(all(grepl("^TCGA_", TCGAcodesAvailable$diseaseCodes)))

# Test no duplicates
expect_equal(
  length(unique(TCGAcodesAvailable$diseaseCodes)),
  nrow(TCGAcodesAvailable)
)

# Test no NA values
expect_false(any(is.na(TCGAcodesAvailable$diseaseCodes)))
expect_false(any(is.na(TCGAcodesAvailable$slide_level_available)))
expect_false(any(is.na(TCGAcodesAvailable$tile_level_available)))
expect_false(any(is.na(TCGAcodesAvailable$hover_available)))

# Test specific known values
expect_true("TCGA_OV" %in% TCGAcodesAvailable$diseaseCodes)

# Test TCGA_OV has HoverNet available (only one with hover_available = TRUE)
ov_row <- TCGAcodesAvailable[TCGAcodesAvailable$diseaseCodes == "TCGA_OV", ]
expect_true(ov_row$hover_available)
expect_equal(sum(TCGAcodesAvailable$hover_available), 1)

# Test all have tile_level available
expect_true(all(TCGAcodesAvailable$tile_level_available))

# Test some known slide_level unavailable
expect_false(
  TCGAcodesAvailable[TCGAcodesAvailable$diseaseCodes == "TCGA_GBM", 
                     "slide_level_available"]
)
expect_false(
  TCGAcodesAvailable[TCGAcodesAvailable$diseaseCodes == "TCGA_KIRC", 
                     "slide_level_available"]
)
