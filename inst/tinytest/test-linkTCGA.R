# Test linkTCGA and helper functions

# Test .rm_tensor_txt
expect_equal(
    ImageFeatureTCGA:::.rm_tensor_txt("tensor(123.)"),
    123
)
expect_equal(
    ImageFeatureTCGA:::.rm_tensor_txt("tensor(456.789.)"),
    456.789
)
expect_equal(
    ImageFeatureTCGA:::.rm_tensor_txt(c("tensor(1.)", "tensor(2.)")),
    c(1, 2)
)

# Test .slide_to_sampleId
expect_equal(
    ImageFeatureTCGA:::.slide_to_sampleId(
        "TCGA-AA-3518-01A-01-BS1.9437f2c5-9f15-4b8d-b95c-01a1bd09b8bd"
    ),
    "TCGA-AA-3518-01A-01-BS1"
)
expect_equal(
    ImageFeatureTCGA:::.slide_to_sampleId(c(
        "TCGA-AA-3518-01A-01-BS1.abc123",
        "TCGA-BB-1234-01A-02-TS2.def456"
    )),
    c("TCGA-AA-3518-01A-01-BS1", "TCGA-BB-1234-01A-02-TS2")
)

# Test .tile_prep_meta
mock_meta <- data.frame(
    slide_name = "slide1",
    tile_id = "tile_1",
    tile_name = "tile_name_1",
    tile_x = "tensor(100.)",
    tile_y = "tensor(200.)"
)
result_meta <- ImageFeatureTCGA:::.tile_prep_meta(mock_meta)
expect_equal(result_meta$tile_x, 100)
expect_equal(result_meta$tile_y, 200)
expect_inherits(result_meta$tile_x, "numeric")
expect_inherits(result_meta$tile_y, "numeric")

# Test .slide_df_to_se
mock_slide_df <- data.frame(
    slideName = c(
        "TCGA-AA-3518-01A-01-BS1.uuid1",
        "TCGA-BB-1234-01A-02-TS2.uuid2"
    ),
    tumorType = c("TCGA_COAD", "TCGA_COAD"),
    fileName = c("file1.csv", "file2.csv"),
    embed_1 = c(0.1, 0.2),
    embed_2 = c(0.3, 0.4),
    embed_3 = c(0.5, 0.6)
)

slide_se <- ImageFeatureTCGA:::.slide_df_to_se(mock_slide_df)

expect_inherits(slide_se, "SummarizedExperiment")
expect_equal(ncol(slide_se), 2)
expect_equal(nrow(slide_se), 3)
expect_true("embeddings" %in% names(SummarizedExperiment::assays(slide_se)))
expect_true("patientIds" %in% names(S4Vectors::metadata(slide_se)))
expect_true("sampleIds" %in% names(S4Vectors::metadata(slide_se)))

# Test .PROV_ORDER and .HOV_ORDER constants
expect_equal(ImageFeatureTCGA:::.PROV_ORDER, c("pipeline", "level", "filename"))
expect_equal(ImageFeatureTCGA:::.HOV_ORDER, c("pipeline", "format", "filename"))

