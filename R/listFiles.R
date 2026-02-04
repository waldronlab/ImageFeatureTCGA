.BASE_URL <- "https://store.cancerdatasci.org"
.PROV_BASE_URL <- paste0(.BASE_URL, "/provgigapath")

#' @name listFiles
#'
#' @title List available HoVerNet and Prov-Giga-Path data for TCGA cancers
#'
#' @description Functions to list available HoverNet and ProvGiga data for TCGA
#'   cancers. HoverNet data is only available for TCGA-OV, while ProvGiga data
#'   is available for multiple TCGA cancer types at slide and tile levels. See
#'   the `TCGAcodesAvailable` dataset for a summary of available data. These
#'   functions return a `data.frame` with filenames and file sizes.
#'
#' @param format `character(1L)` One of "geojson", "h5ad", "json", or "thumb"
#'   specifying the desired HoverNet data format. Default is "h5ad".
#'
#' @param level `character(1L)` One of "slide_level" or "tile_level" specifying
#'   the desired ProvGiga data level. Default is "slide_level".
#'
#' @returns `listHoverNet`,`listProvGiga`: A `tibble` listing available HoverNet
#'   or ProvGigaPath files with `Filename`, `Modified`, and `Size` columns.
#'
#' @examplesIf interactive()
#' ## List available HoverNet data for TCGA-OV
#' listHoverNet(format = "h5ad")
#' @export
listHoverNet <- function(
    format = c("geojson", "h5ad", "json", "thumb")
) {
    format <- match.arg(format)
    hovernet_url <-
        paste(.BASE_URL, "hovernet", format, "", sep = "/")
    table <- .see_more_table(hovernet_url)
    table[!grepl("^\\.\\.", table[["Filename"]]), ]
}

#' @rdname listFiles
#'
#' @importFrom BiocBaseUtils isScalarCharacter
#'
#' @examplesIf interactive()
#' ## List available ProvGiga slide-level data for TCGA-BRCA
#' listProvGiga(level = "slide_level")
#' @export
listProvGiga <- function(
    level = c("slide_level", "tile_level")
) {
    level <- match.arg(level)

    tumor_type_url <- paste(
        .PROV_BASE_URL, level, "", sep = "/"
    )
    table <- .see_more_table(tumor_type_url)
    table[!grepl("^\\.\\.", table[["Filename"]]), ]
}

.CATALOG_COL_TYPES <- "ccccccccccccccccccccccddc"

#' @rdname listFiles
#'
#' @description The `getCatalog` function retrieves a catalog of all available
#'   HoVerNet and ProvGigaPath files, including filenames, sizes, pipelines
#'   used, tumor types, and data levels.
#'
#' @param pipeline `character()` One or both "hovernet" and/or "provgigapath"
#'   specifying which pipeline(s) to include in the catalog. Default includes
#'   both.
#'
#' @param format `character()` One or more of "csv", "thumb", "h5ad", "geojson",
#'   or "json" specifying which file formats to include in the catalog. Default
#'   includes all.
#'
#' @param redownload `logical(1L)` Whether to redownload the catalog file even
#'   if it is already cached locally. Default is `FALSE`.
#'
#' @returns `getCatalog`: A `tibble` containing the full catalog of available
#'   files for the specified pipeline(s).
#'
#' @examplesIf interactive()
#' ## Get the full catalog of available files
#' getCatalog(pipeline = c("hovernet", "provgigapath"), format = "h5ad")
#' @export
getCatalog <-
    function(
        pipeline = c("hovernet", "provgigapath"),
        format = c("csv", "thumb", "h5ad", "geojson", "json"),
        redownload = FALSE
    )
{
    pipeline <- match.arg(pipeline, several.ok = TRUE)
    format <- match.arg(format, several.ok = TRUE)
    catalog <- .download_catalog(redownload = redownload) |>
        readr::read_tsv(col_types = .CATALOG_COL_TYPES)
    in_pipe <- catalog[["pipeline"]] %in% pipeline
    in_format <- catalog[["format"]] %in% format
    catalog[in_pipe & in_format, ]
}

#' @importFrom httr2 request req_headers req_perform resp_body_json
.download_catalog <- function(redownload) {
    resp <- request("https://zenodo.org/api/records/17981132") |>
        req_headers(
            Accept = "application/json"
        ) |>
        req_perform() |>
        resp_body_json()

    cache <- getOption(
        "BiocFileCache.cache", BiocFileCache::getBFCOption("CACHE")
    )
    bfc <- BiocFileCache::BiocFileCache(cache = cache)
    .cache_url_file(
        resp$files[[1L]]$links$self, redownload = redownload, bfc = bfc
    )
}

#' @rdname listFiles
#'
#' @param catalog A `tibble` as returned by `getCatalog()`.
#'
#' @returns `getFileURLs`: A `character()` vector of full URLs for the files
#'   listed in the provided catalog.
#'
#' @examplesIf interactive()
#' ## Get file URLs from the catalog
#' getCatalog(pipeline = "hovernet", format = "h5ad") |>
#'     dplyr::slice(1:10) |>
#'     getFileURLs()
#' @export
getFileURLs <- function(catalog) {
    paste(.BASE_URL, catalog[["fullpath"]], sep = "/")
}
