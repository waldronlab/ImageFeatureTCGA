#' @importFrom BiocBaseUtils isScalarCharacter
.is_url <- function(url) {
    stopifnot(
        isScalarCharacter(url)
    )
    grepl("^https?://|^ftp://", url)
}

.cache_url_file <- function(url, redownload = FALSE) {
    cache <- getOption(
        "BiocFileCache.cache", BiocFileCache::getBFCOption("CACHE")
    )
    bfc <- BiocFileCache::BiocFileCache(cache = cache)
    bquery <- BiocFileCache::bfcquery(bfc, url, "rname", exact = TRUE)
    cached <- identical(nrow(bquery), 1L)

    if (!redownload && cached)
        return(bquery[["rpath"]])

    part_url <- gsub(paste0(.BASE_URL, "/"), "", url)
    destfile <- file.path(cache, part_url)
    destfolder <- dirname(destfile)
    if (!dir.exists(destfolder))
        dir.create(destfolder, recursive = TRUE, showWarnings = FALSE)
    file <- curl::curl_download(url = url, destfile = destfile)
    if (!cached)
        BiocFileCache::bfcadd(
            x = bfc,
            rname = url,
            fpath = file,
            rtype = "local",
            action = "asis",
            fname = "exact",
            exact = TRUE
        )
    else
        file
}

.url_query <- function(bfc, urls) {
    lapply(
        urls,
        function(url) {
            BiocFileCache::bfcquery(bfc, url, "rname", exact = TRUE)
        }
    )
}

.is_cached <- function(qframe) {
    vapply(qframe, nrow, integer(1L)) == 1L
}

.rpath_cache <- function(qframe) {
    vapply(qframe, `[[`, character(1L), "rpath")
}

.cache_url_files <- function(urls, redownload = FALSE, parallel = TRUE) {
    checkInstalled("curl")
    checkInstalled("BiocFileCache")
    if (parallel) {
        cache <- getOption(
            "BiocFileCache.cache", BiocFileCache::getBFCOption("CACHE")
        )
        bfc <- BiocFileCache::BiocFileCache(cache = cache)
        queries <- .url_query(bfc, urls)
        cached <- .is_cached(queries)
        locals <- vector("list", length(urls))

        if (!redownload)
            locals[cached] <- .rpath_cache(queries[cached])
        urls <- urls[!cached | redownload]
        if (length(urls)) {
            part_urls <- gsub(paste0(.BASE_URL, "/"), "", urls)
            destfiles <- file.path(
                BiocFileCache::getBFCOption("CACHE"), part_urls
            )
            destfolders <- dirname(destfiles) |>
                unique()
            dexist <- dir.exists(destfolders)
            if (!all(dexist))
                vapply(
                    destfolders[!dexist],
                    dir.create,
                    logical(1L),
                    recursive = TRUE,
                    showWarnings = FALSE
                )
            output <- curl::multi_download(
                urls = urls,
                destfiles = destfiles
            )
            successframe <- output[output[["success"]], , drop = FALSE]
            successurls <- urls[output[["success"]]]
            successfiles <- successframe[["destfile"]]

            locals <- BiocParallel::bpmapply(
                function(bfc, url, file, cached, bfcid) {
                    if (!cached)
                        BiocFileCache::bfcadd(
                            x = bfc,
                            rname = url,
                            fpath = file,
                            rtype = "local",
                            action = "asis",
                            fname = "exact",
                            exact = TRUE
                        )
                    else
                        file
                },
                url = successurls,
                file = successfiles,
                cached = cached,
                MoreArgs = list(bfc = bfc),
                SIMPLIFY = FALSE
            )
        }
        unlist(locals)
    } else {
        vapply(
            urls,
            function(url) {
                .cache_url_file(
                    url = url,
                    redownload = redownload
                )
            },
            character(1L)
        )
    }
}

#' @importFrom rvest html_nodes html_table html_element html_attr read_html
#' @importFrom httr2 request req_perform resp_body_string
.see_more_table <- function(u24_url, verbose = TRUE) {
    results <- list()
    current_url <- u24_url
    page_count <- 1

    while (!is.null(current_url)) {
        page_html <- request(current_url) |>
            req_perform() |>
            resp_body_string() |>
            read_html()

        table_data <- page_html |>
            html_nodes("table") |>
            html_table(fill = TRUE)

        results[[page_count]] <- table_data

        cursor_node <- page_html |>
            html_element("a:contains('see more')")

        if (!is.na(cursor_node)) {
            query_string <- html_attr(cursor_node, "href")
            current_url <- paste0(u24_url, query_string)
            page_count <- page_count + 1
        } else {
            current_url <- NULL
            if (verbose)
                message("Total pages fetched: ", page_count)
        }
    }
    dplyr::bind_rows(
        unlist(results, recursive = FALSE)
    )
}

.import_slide_level <- function(
    prov_path, tumorType, fileName, layer = "last_layer_embed", ...
) {
    df <- readr::read_csv(prov_path, show_col_types = FALSE)
    embedding <- df[[layer]][1L] |>
        gsub("tensor\\(\\[\\[|\\]\\]\\)", "", x = _) |>
        gsub("\\n", "", x = _) |>
        read.table(text = _, sep = ",")

    tibble::tibble(
        slideName = df[["slide_name"]],
        tumorType = tumorType,
        fileName = fileName,
        embedding
    )
}

.import_tile_level <- function(prov_path, tumorType, fileName, ...) {
    df <- readr::read_csv(prov_path, show_col_types = FALSE)
    tibble::tibble(
        df,
        tumorType = tumorType,
        fileName = fileName
    )
}

.extract_tcgabcode <- function(file_path) {
    stopifnot(
        isScalarCharacter(file_path)
    )
    utils::head(
        strsplit(basename(file_path), "\\.")[[1L]],
        1L
    )
}


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
#' @importFrom data.table data.table as.data.table setnames .N := .I
#' @importFrom dplyr mutate
#' @importFrom methods is
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

    # Create xmin, xmax, ymin, ymax columns using column names
    dt_tiles$xmin <- dt_tiles[[tile_x]]
    dt_tiles$xmax <- dt_tiles[[tile_x]] + tile_size
    dt_tiles$ymin <- dt_tiles[[tile_y]]
    dt_tiles$ymax <- dt_tiles[[tile_y]] + tile_size

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
    # Use .SD to work with the tile_id column name
    tile_counts <- assigned[, .N, by = c(tile_id, "cell_type")]
    data.table::setnames(tile_counts, "cell_type", "cell_type_label")

    # Select the most frequent cell type in each tile
    tile_dominant <- tile_counts[
        tile_counts[, .I[which.max(N)], by = tile_id]$V1
    ]
    data.table::setnames(tile_dominant, "cell_type_label", "dominant_cell_type")

    # Step 6: Merge cell type info back to tiles
    tiles_with_nuclei <- merge(
        as.data.frame(tiles),
        as.data.frame(tile_counts),
        by = tile_id,
        all.x = TRUE
    )

    tiles_with_nuclei <- merge(
        tiles_with_nuclei,
        as.data.frame(tile_dominant)[, c(tile_id, "dominant_cell_type")],
        by = tile_id,
        all.x = TRUE
    )

    # Return results
    list(
        tiles_with_nuclei = tiles_with_nuclei,
        tiles_dominant = as.data.frame(tile_dominant),
        scale_factor = sf_list
    )
}
