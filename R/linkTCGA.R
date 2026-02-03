.PROV_ORDER <- c("pipeline", "level", "filename")
.HOV_ORDER <- c("pipeline", "format", "filename")

#' Link MultiAssayExperiment object to TCGA data
#'
#' @param MultiAssayExperiment A `MultiAssayExperiment` object containing sample
#'   metadata with TCGA barcodes.
#'
#' @importFrom S4Vectors DataFrame
#'
#' @examplesIf interactive()
#' library(curatedTCGAData)
#' coad <- curatedTCGAData(
#'     diseaseCode = "COAD",
#'     assays = "RNASeq2GeneNorm",
#'     version = "2.1.1",
#'     dry.run = FALSE
#' )
#' coad_sub <- coad[, 1:3, ]
#' catalog <- getCatalog(pipeline = "provgigapath", format = "csv")
#' linkTCGA(coad_sub, catalog, parallel = FALSE)
#' @export
linkTCGA <- function(
    MultiAssayExperiment, catalog, redownload = FALSE, parallel = TRUE
) {
    tcgabcodes <- TCGAutils::TCGAbarcode(catalog[["tcga_barcode"]])
    catalog <-
        catalog[tcgabcodes %in% rownames(colData(MultiAssayExperiment)), ]
    catalog[["url"]] <- getFileURLs(catalog)
    resdata <- ProvGigaList(
        catalog[["url"]],
        is_url = TRUE,
        levels = catalog[["level"]],
        parallel = parallel
    ) |>
        import(redownload = redownload, parallel = parallel)
    slide_assay <- .slide_df_to_se(resdata[["slide_level"]])
    slide_sm <- DataFrame(
        assay = "slide_assay",
        primary = metadata(slide_assay)[["patientIds"]],
        colname = metadata(slide_assay)[["sampleIds"]]
    )
    tile_assay <- .tile_df_to_bumpy_se(resdata[["tile_level"]])
    tile_sm <- DataFrame(
        assay = "tile_assay",
        primary = metadata(tile_assay)[["patientIds"]],
        colname = metadata(tile_assay)[["sampleIds"]]
    )
    sampmap <- rbind(slide_sm, tile_sm)
    c(
        MultiAssayExperiment,
        slide_assay = slide_assay,
        tile_assay = tile_assay,
        sampleMap = sampmap
    )
}

#' @importFrom SummarizedExperiment SummarizedExperiment
.slide_df_to_se <- function(sdf) {
    sampleIds <- .slide_to_sampleId(sdf[["slideName"]])
    patientIds <- TCGAutils::TCGAbarcode(sampleIds)
    metadata <- append(
        as.list(sdf[, c("slideName", "tumorType", "fileName")]),
        list(
            patientIds = patientIds,
            sampleIds = sampleIds
        )
    )
    embeddings <-
        sdf[-which(names(sdf) %in% c("slideName", "tumorType", "fileName"))] |>
        as.matrix() |>
        t()
    dimnames(embeddings) <- list(
        NULL,
        sampleIds
    )
    SummarizedExperiment(
        assays = list(embeddings = embeddings),
        metadata = metadata
    )
}

.tile_df_to_bumpy_se <- function(tdf) {
    meta <- c("slide_name", "tile_id", "tile_name", "tile_x", "tile_y")
    mdat <- .tile_prep_meta(tdf[, meta])
    sampleIds <- .slide_to_sampleId(tdf[["slide_name"]])
    patientIds <- TCGAutils::TCGAbarcode(sampleIds)
    metadata <- list(
        metadata = mdat, sampleIds = sampleIds, patientIds = patientIds
    )
    ncols <- ncol(tdf) - length(meta)
    nsamps <- length(unique(sampleIds))
    enames <- tdf[["slide_name"]] |>
        as.factor() |>
        levels() |>
        .slide_to_sampleId()
    split(
        tdf[, -which(colnames(tdf) %in% meta)],
        tdf[["slide_name"]]
    ) |>
        lapply(function(x) unname(x) |> IRanges::NumericList()) |>
        unname() |>
        do.call(c, args = _) |>
        BumpyMatrix::BumpyMatrix(
            dim = c(ncols, length(enames)),
            dimnames = list(NULL, enames)
        ) |>
        list(tiles = _) |>
        SummarizedExperiment(
            metadata = metadata
        )
}

.tile_prep_meta <- function(metadf) {
    metadf[["tile_x"]] <- .rm_tensor_txt(metadf[["tile_x"]])
    metadf[["tile_y"]] <- .rm_tensor_txt(metadf[["tile_y"]])
    metadf
}

.rm_tensor_txt <- function(txt) {
    gsub("tensor\\((.*)\\.\\)", "\\1", x = txt) |>
        as.numeric()
}

.slide_to_sampleId <- function(txt) {
    vapply(strsplit(txt, "\\."), `[[`, character(1L), 1L)
}
