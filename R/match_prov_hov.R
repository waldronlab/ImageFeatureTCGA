#' Match HoverNet Nuclei to ProvGigaPath Tiles
#'
#' @description Assigns HoverNet nuclei to ProvGigaPath tiles by computing a
#'   scale factor to align coordinate systems, then performing spatial matching.
#'   Returns tile-level cell type counts and dominant cell types.
#'
#' @param hovernet A `SpatialExperiment` or `SpatialFeatureExperiment` object
#'   containing HoverNet nuclei data with spatial coordinates and cell type
#'   information.
#' @param tiles A `data.frame` or `tibble` containing ProvGigaPath tile data
#'   with tile coordinates and embeddings.
#' @param tile_size Numeric. Size of tiles in pixels. Default is `224`.
#' @param cell_x Character. Name of the x-coordinate column in HoverNet data.
#'   Default is `"x_centroid"` for h5ad data.
#' @param cell_y Character. Name of the y-coordinate column in HoverNet data.
#'   Default is `"y_centroid"` for h5ad data.
#' @param tile_x Character. Name of the tile x-coordinate column in tiles data.
#'   Default is `"tile_x"`.
#' @param tile_y Character. Name of the tile y-coordinate column in tiles data.
#'   Default is `"tile_y"`.
#' @param tile_id Character. Name of the tile ID column in tiles data.
#'   Default is `"tile_id"`.
#' @param cell_type Character. Name of the cell type column in HoverNet data.
#'   Default is `"type"`.
#'
#' @return A list with three elements:
#'   \describe{
#'     \item{tiles_with_nuclei}{A `data.frame` with original tile data plus
#'       cell type counts (`N`) and labels (`cell_type_label`) for each tile.}
#'     \item{tiles_dominant}{A `data.frame` with tile IDs and their dominant
#'       (most frequent) cell type.}
#'     \item{scale_factor}{A list containing the computed scale factor, and
#'       separate x and y scale factors.}
#'   }
#'
#' @details The function performs the following steps:
#'   1. Computes a scale factor to align nuclei coordinates with tile coordinates
#'   2. Scales nuclei coordinates using the computed scale factor
#'   3. Creates bounding boxes for each tile based on tile_size
#'   4. Assigns nuclei to tiles using spatial overlap
#'   5. Counts cell types per tile
#'   6. Identifies the dominant cell type for each tile
#'   7. Merges results back to the original tiles data
#'
#'   The scale factor is computed as the mean of x and y scale factors, where
#'   each is the ratio of coordinate ranges between nuclei and tiles.
#'
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom data.table data.table as.data.table setnames .N :=
#' @importFrom dplyr mutate
#' @importFrom methods is
#' @importFrom grDevices rgb
#'
#' @examples
#' \dontrun{
#' # Basic usage with h5ad HoverNet data
#' result <- matchHoverNetToTiles(hn_spe, pca_tiles)
#'
#' # Access results
#' tiles_matched <- result$tiles_with_nuclei
#' dominant_types <- result$tiles_dominant
#' scale_info <- result$scale_factor
#'
#' # Custom parameters
#' result <- matchHoverNetToTiles(
#'   hn_spe,
#'   pca_tiles,
#'   tile_size = 256,
#'   tile_x = "x_coord",
#'   tile_y = "y_coord"
#' )
#' }
#'
#' @export
matchHoverNetToTiles <- function(
    hovernet,
    tiles,
    tile_size = 224,
    cell_x = "x_centroid",
    cell_y = "y_centroid",
    tile_x = "tile_x",
    tile_y = "tile_y",
    tile_id = "tile_id",
    cell_type = "type"
) {
    
    x <- y <- xmin <- xmax <- ymin <- ymax <- x1 <- y1 <- N <- NULL
    # Validate input
    if (!is(hovernet, "SpatialExperiment") && 
        !is(hovernet, "SpatialFeatureExperiment")) {
        stop("'hovernet' must be a SpatialExperiment or ",
             "SpatialFeatureExperiment object.")
    }
    
    if (!is.data.frame(tiles)) {
        stop("'tiles' must be a data.frame or tibble.")
    }
    
    # Extract nuclei metadata
    coords <- spatialCoords(hovernet)
    cell_data <- colData(hovernet)
    
    # Check required columns
    if (!all(c(cell_x, cell_y) %in% colnames(coords))) {
        stop("Spatial coordinates '", cell_x, "' and '", cell_y, 
             "' not found in spatialCoords.")
    }
    
    if (!cell_type %in% colnames(cell_data)) {
        stop("Cell type column '", cell_type, "' not found in colData.")
    }
    
    if (!all(c(tile_x, tile_y, tile_id) %in% colnames(tiles))) {
        stop("Required tile columns not found: ", 
             paste(c(tile_x, tile_y, tile_id), collapse = ", "))
    }
    
    # Create cell metadata data.frame
    cell_meta <- data.frame(
        x = coords[, cell_x],
        y = coords[, cell_y],
        type = cell_data[[cell_type]]
    )
    
    # Step 1: Compute scale factor
    sx <- diff(range(cell_meta$x, na.rm = TRUE)) /
        diff(range(tiles[[tile_x]], na.rm = TRUE))
    sy <- diff(range(cell_meta$y, na.rm = TRUE)) /
        diff(range(tiles[[tile_y]], na.rm = TRUE))
    scale_factor <- mean(c(sx, sy), na.rm = TRUE)
    
    sf_list <- list(
        scale_factor = scale_factor,
        sx = sx,
        sy = sy
    )
    
    # Step 2: Scale nuclei coordinates
    cell_meta_scaled <- cell_meta |>
        dplyr::mutate(
            x = x / scale_factor,
            y = y / scale_factor
        )
    
    # Step 3: Prepare tiles and nuclei data.tables
    dt_tiles <- data.table::as.data.table(tiles)
    dt_tiles[, `:=`(
        xmin = get(tile_x),
        xmax = get(tile_x) + tile_size,
        ymin = get(tile_y),
        ymax = get(tile_y) + tile_size
    )]
    
    # Nuclei table: each nucleus as a point
    dt_nuc <- data.table::data.table(
        x1 = cell_meta_scaled$x,
        x2 = cell_meta_scaled$x,
        y1 = cell_meta_scaled$y,
        y2 = cell_meta_scaled$y,
        cell_type = cell_meta_scaled$type
    )
    
    # Step 4: Assign nuclei to tiles using non-equi join
    assigned <- dt_tiles[
        dt_nuc,
        on = .(xmin <= x1, xmax >= x1, ymin <= y1, ymax >= y1),
        nomatch = 0L
    ]
    
    # Step 5: Compute per-tile cell type counts
    tile_counts <- assigned[, .N, by = c(tile_id, "cell_type")]
    data.table::setnames(tile_counts, "cell_type", "cell_type_label")
    
    # Select the most frequent cell type in each tile
    tile_dominant <- tile_counts[
        tile_counts[, .I[which.max(N)], by = get(tile_id)]$V1
    ]
    data.table::setnames(tile_dominant, "cell_type_label", "dominant_cell_type")
    
    # Step 6: Merge cell type info back to tiles
    tiles_with_nuclei <- merge(
        tiles,
        tile_counts,
        by = tile_id,
        all.x = TRUE
    )
    
    tiles_with_nuclei <- merge(
        tiles_with_nuclei,
        tile_dominant[, c(tile_id, "dominant_cell_type"), with = FALSE],
        by = tile_id,
        all.x = TRUE
    )
    
    # Return results
    list(
        tiles_with_nuclei = as.data.frame(tiles_with_nuclei),
        tiles_dominant = as.data.frame(tile_dominant),
        scale_factor = sf_list
    )
}
