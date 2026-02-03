# Test matchHoverNetToTiles

# Create mock SpatialExperiment for HoverNet data
library(SpatialExperiment)
library(S4Vectors)

# Mock nuclei data (100 nuclei)
set.seed(123)
n_nuclei <- 100
coords_mat <- matrix(
    c(runif(n_nuclei, 0, 10000), runif(n_nuclei, 0, 10000)),
    ncol = 2,
    dimnames = list(NULL, c("x_centroid", "y_centroid"))
)

cell_types <- sample(0:5, n_nuclei, replace = TRUE)

# Create minimal SpatialExperiment
mock_spe <- SpatialExperiment(
    assays = list(counts = matrix(0, nrow = 1, ncol = n_nuclei)),
    spatialCoords = coords_mat,
    colData = DataFrame(type = cell_types)
)

# Mock tiles data (25 tiles in 5x5 grid)
tile_coords <- expand.grid(
    tile_x = seq(0, 800, by = 224),
    tile_y = seq(0, 800, by = 224)
)
mock_tiles <- data.frame(
    tile_id = paste0("tile_", seq_len(nrow(tile_coords))),
    tile_x = tile_coords$tile_x,
    tile_y = tile_coords$tile_y,
    embedding_1 = runif(nrow(tile_coords))
)

# Test basic functionality
result <- matchHoverNetToTiles(mock_spe, mock_tiles)

expect_inherits(result, "list")
expect_equal(length(result), 3)
expect_true(all(c("tiles_with_nuclei", "tiles_dominant", "scale_factor") %in% names(result)))

# Test tiles_with_nuclei output
expect_inherits(result$tiles_with_nuclei, "data.frame")
expect_true("cell_type_label" %in% names(result$tiles_with_nuclei))
expect_true("N" %in% names(result$tiles_with_nuclei))
expect_true("dominant_cell_type" %in% names(result$tiles_with_nuclei))
expect_true("tile_id" %in% names(result$tiles_with_nuclei))

# Test tiles_dominant output
expect_inherits(result$tiles_dominant, "data.frame")
expect_true("dominant_cell_type" %in% names(result$tiles_dominant))
expect_true("N" %in% names(result$tiles_dominant))

# Test scale_factor output
expect_inherits(result$scale_factor, "list")
expect_true(all(c("scale_factor", "sx", "sy") %in% names(result$scale_factor)))
expect_inherits(result$scale_factor$scale_factor, "numeric")
expect_true(result$scale_factor$scale_factor > 0)

# Test custom tile_size
result_256 <- matchHoverNetToTiles(mock_spe, mock_tiles, tile_size = 256)
expect_inherits(result_256, "list")

# Test input validation - wrong hovernet class
expect_error(
    matchHoverNetToTiles(data.frame(), mock_tiles),
    "must be a SpatialExperiment"
)

# Test input validation - wrong tiles class
expect_error(
    matchHoverNetToTiles(mock_spe, "not_a_dataframe"),
    "must be a data.frame"
)

# Test input validation - missing tile columns
bad_tiles <- data.frame(x = 1, y = 2)
expect_error(
    matchHoverNetToTiles(mock_spe, bad_tiles),
    "Required tile columns not found"
)

# Test input validation - missing cell type column
mock_spe_no_type <- mock_spe
colData(mock_spe_no_type)$type <- NULL
expect_error(
    matchHoverNetToTiles(mock_spe_no_type, mock_tiles),
    "not found in colData"
)

# Test custom column names
mock_tiles_custom <- mock_tiles
names(mock_tiles_custom) <- c("id", "x", "y", "emb")
result_custom <- matchHoverNetToTiles(
    mock_spe, 
    mock_tiles_custom,
    tile_x = "x",
    tile_y = "y",
    tile_id = "id"
)
expect_inherits(result_custom, "list")

