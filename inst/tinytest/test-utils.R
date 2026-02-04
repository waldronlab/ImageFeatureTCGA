# Test utils.R functions

# Test .is_url
expect_true(ImageFeatureTCGA:::.is_url("https://example.com"))
expect_true(ImageFeatureTCGA:::.is_url("http://example.com"))
expect_true(ImageFeatureTCGA:::.is_url("ftp://example.com"))
expect_false(ImageFeatureTCGA:::.is_url("/path/to/file.csv"))
expect_false(ImageFeatureTCGA:::.is_url("file.csv"))

# Test .is_url error on non-scalar
expect_error(
    ImageFeatureTCGA:::.is_url(c("https://a.com", "https://b.com"))
)
expect_error(
    ImageFeatureTCGA:::.is_url(123)
)

# Test .extract_tcgabcode
expect_equal(
    ImageFeatureTCGA:::.extract_tcgabcode(
        "/path/to/TCGA-AA-3518-01A-01-BS1.9437f2c5-9f15-4b8d-b95c-01a1bd09b8bd.csv"
    ),
    "TCGA-AA-3518-01A-01-BS1"
)
expect_equal(
    ImageFeatureTCGA:::.extract_tcgabcode("TCGA-BB-1234.uuid.csv.gz"),
    "TCGA-BB-1234"
)

# Test .extract_tcgabcode error on non-scalar
expect_error(
    ImageFeatureTCGA:::.extract_tcgabcode(c("file1.csv", "file2.csv"))
)

# Test .is_cached
mock_query_cached <- list(
    data.frame(rpath = "/path/to/file1"),
    data.frame(rpath = "/path/to/file2")
)
mock_query_not_cached <- list(
    data.frame(rpath = character(0)),
    data.frame(rpath = "/path/to/file2")
)
expect_equal(
    ImageFeatureTCGA:::.is_cached(mock_query_cached),
    c(TRUE, TRUE)
)
expect_equal(
    ImageFeatureTCGA:::.is_cached(mock_query_not_cached),
    c(FALSE, TRUE)
)

# Test .rpath_cache
expect_equal(
    ImageFeatureTCGA:::.rpath_cache(mock_query_cached),
    c("/path/to/file1", "/path/to/file2")
)

# Test .import_slide_level (requires temp file)
mock_slide_csv <- tempfile(fileext = ".csv")
writeLines(
    c(
        "slide_name,last_layer_embed",
        "TCGA-AA-3518-01A.uuid,\"tensor([[0.1, 0.2, 0.3]])\""
    ),
    mock_slide_csv
)
slide_result <- ImageFeatureTCGA:::.import_slide_level(
    prov_path = mock_slide_csv,
    tumorType = "TCGA_COAD",
    fileName = "test.csv"
)
expect_inherits(slide_result, "data.frame")
expect_true("slideName" %in% names(slide_result))
expect_true("tumorType" %in% names(slide_result))
expect_true("fileName" %in% names(slide_result))
expect_equal(slide_result$tumorType, "TCGA_COAD")
expect_equal(ncol(slide_result), 6)  
unlink(mock_slide_csv)

# Test .import_tile_level (requires temp file)
mock_tile_csv <- tempfile(fileext = ".csv")
writeLines(
    c(
        "slide_name,tile_id,tile_x,tile_y,embed_1",
        "slide1,tile_1,0,0,0.5",
        "slide1,tile_2,224,0,0.6"
    ),
    mock_tile_csv
)
tile_result <- ImageFeatureTCGA:::.import_tile_level(
    prov_path = mock_tile_csv,
    tumorType = "TCGA_COAD",
    fileName = "test_tile.csv"
)
expect_inherits(tile_result, "data.frame")
expect_true("tumorType" %in% names(tile_result))
expect_true("fileName" %in% names(tile_result))
expect_equal(nrow(tile_result), 2)
expect_equal(tile_result$tumorType[1], "TCGA_COAD")
unlink(mock_tile_csv)

# Test .PROV_ORDER and .HOV_ORDER constants
expect_equal(
    ImageFeatureTCGA:::.PROV_ORDER, 
    c("pipeline", "level", "filename")
)
expect_equal(
    ImageFeatureTCGA:::.HOV_ORDER, 
    c("pipeline", "format", "filename")
)


