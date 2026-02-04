#' Plot HoverNet Segmentation Overlay with Thumbnail
#'
#' @description Creates a side-by-side visualization of the original tissue
#'   thumbnail image and the HoverNet cell segmentation with colored cell type
#'   labels. The function automatically retrieves the thumbnail image associated
#'   with a HoverNet JSON file and overlays the segmentation data.
#'
#' @param hovernet A `SpatialExperiment` or `SpatialFeatureExperiment` object
#'   created by importing a HoverNet JSON file, or a path/URL to a HoverNet
#'   JSON file.
#' @param json_path Optional. Path or URL to the HoverNet JSON file. Only
#'   required if `hovernet` is a `SpatialExperiment` object and the original
#'   JSON path is not stored in metadata. Default is `NULL`.
#' @param title Optional. Title for the combined plot. If `NULL`, uses the
#'   basename of the JSON file. Default is `NULL`.
#' @param point_size Numeric value for the size of points in the segmentation
#'   plot. Default is `0.01`.
#' @param legend_point_size Numeric value for the size of points in the legend.
#'   Default is `2`.
#' @param color_palette Optional. A named vector of colors for cell types. If
#'   `NULL`, uses the colors from the HoverNet type map or RColorBrewer.
#'   Default is `NULL`.
#' @param ncol Number of columns for the side-by-side plots. Default is `2`.
#' @param rel_widths Relative widths of the image and segmentation plots.
#'   Default is `c(1, 1.08)`.
#' @param title_size Font size for the main title. Default is `25`.
#' @param flip_image Logical. Whether to flip the image vertically to match
#'   coordinate system. Default is `TRUE`.
#'
#' @return A `ggplot` object combining the thumbnail image and segmentation
#'   overlay with an optional title.
#'
#' @details The function performs the following steps:
#'   1. If `hovernet` is a file path, imports it using `HoverNet()`
#'    and `import()`
#'   2. Retrieves the associated thumbnail PNG using `importHoverNetThumbnail()`
#'   3. Creates a segmentation plot colored by cell type
#'   4. Combines the thumbnail and segmentation side-by-side
#'   5. Adds an optional title
#'
#'    The thumbnail is automatically downloaded and cached if a URL is provided.
#'    The segmentation uses the `type` column from the `colData` and colors are
#'    taken from the `type_map` metadata if available.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual guides
#'   guide_legend theme_void theme coord_fixed element_blank
#' @importFrom cowplot plot_grid ggdraw draw_image draw_label
#' @importFrom BiocBaseUtils checkInstalled
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom S4Vectors metadata
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # From a JSON file path
#' json_file <- "path/to/sample.json.gz"
#' plotHoverNetOverlay(json_file, title = "Sample Tissue")
#'
#' # From a URL
#' json_url <- paste0(
#'   "https://store.cancerdatasci.org/hovernet/TCGA_OV/json/",
#'   "TCGA-VG-A8LO-01A-01-DX1.B39A4D64-82A1-4A04-8AB6-918F3058B83B.json.gz"
#' )
#' plotHoverNetOverlay(json_url)
#'
#' # From a SpatialExperiment object
#' hn_spe <- HoverNet(json_file) |> import()
#' plotHoverNetOverlay(hn_spe, json_path = json_file)
#'
#' # Custom colors
#' custom_colors <- c(
#'   "neopla" = "#FF0000",
#'   "inflam" = "#00FF00",
#'   "connec" = "#0000FF"
#' )
#' plotHoverNetOverlay(json_file, color_palette = custom_colors)
#' }
#'
#' @export
plotHoverNetOverlay <- function(
    hovernet,
    json_path = NULL,
    title = NULL,
    point_size = 0.01,
    legend_point_size = 2,
    color_palette = NULL,
    ncol = 2,
    rel_widths = c(1, 1.08),
    title_size = 25,
    flip_image = TRUE
) {

    checkInstalled("magick")

    # Handle input: if hovernet is a path/URL, import it
    if (is.character(hovernet) && length(hovernet) == 1) {
        json_path <- hovernet
        hovernet <- HoverNet(json_path, outClass = "SpatialExperiment") |>
            import()
    }

    # Validate input
    if (!is(hovernet, "SpatialExperiment") &&
        !is(hovernet, "SpatialFeatureExperiment")) {
        stop("'hovernet' must be a SpatialExperiment, ",
            "SpatialFeatureExperiment, or a path to a HoverNet JSON file.")
    }

    # Get JSON path from metadata if not provided
    if (is.null(json_path)) {
        json_path <- metadata(hovernet)$json_path
        if (is.null(json_path)) {
            stop("'json_path' must be provided when 'hovernet' is a ",
                    "SpatialExperiment object without metadata$json_path.")
        }
    }

    # Import thumbnail
    png_img <- importHoverNetThumbnail(json_path)
    if (is.null(png_img)) {
        stop("Could not retrieve thumbnail image for: ", json_path)
    }

    # Prepare data for plotting
    gg <- data.frame(
        spatialCoords(hovernet),
        colData(hovernet)
    )

    # Determine color palette
    if (is.null(color_palette)) {
        # Try to use colors from type_map in metadata
        type_map <- metadata(hovernet)$type_map
        if (!is.null(type_map) &&
            all(c("label", "R", "G", "B") %in% names(type_map))) {
            color_palette <- rgb(
                type_map$R / 255,
                type_map$G / 255,
                type_map$B / 255
            )
            names(color_palette) <- type_map$label
        } else {
            n_types <- length(unique(gg$label))
            color_palette <- RColorBrewer::brewer.pal(
                min(n_types, 12), "Paired"
            )[seq_len(n_types)]
        }
    }

    # Create segmentation plot
    p_feat <- ggplot2::ggplot(gg, ggplot2::aes(x, y, col = label)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_manual(
            values = color_palette,
            name = "Cell Type"
        ) +
        ggplot2::guides(
            col = ggplot2::guide_legend(
                override.aes = list(size = legend_point_size)
            )
        ) +
        ggplot2::theme_void()

    # Create version without legend for side-by-side comparison
    p_feat_nolegend <- p_feat +
        ggplot2::theme(legend.position = "none") +
        ggplot2::coord_fixed(ratio = 1)

    # Extract legend
    legend <- cowplot::get_legend(
        p_feat +
            ggplot2::theme(
                legend.position = "right",
                legend.box.margin = ggplot2::margin(0, 0, 0, 10)
            )
    )

    # Prepare image
    img <- magick::image_read(png_img)
    if (flip_image) {
        img <- magick::image_flip(img)
    }

    # Create image plot
    p_img <- cowplot::ggdraw() + cowplot::draw_image(img)

    # Combine plots side by side
    combo_plots <- cowplot::plot_grid(
        p_img, p_feat_nolegend,
        ncol = ncol,
        align = "hv",
        axis = "tblr",
        rel_widths = rel_widths
    )

    # Add legend on the right without changing plot proportions
    combo <- cowplot::plot_grid(
        combo_plots, legend,
        ncol = 2,
        rel_widths = c(1, 0.15)
    )

    # Add title if provided
    if (is.null(title)) {
        title <- tools::file_path_sans_ext(
            tools::file_path_sans_ext(basename(json_path))
        )
    }

    if (!is.null(title) && nchar(title) > 0) {
        title_plot <- cowplot::ggdraw() +
            cowplot::draw_label(
                title,
                fontface = "bold",
                x = 0.5,
                hjust = 0.5,
                size = title_size
            )

        final <- cowplot::plot_grid(
            title_plot, combo,
            ncol = 1,
            rel_heights = c(0.08, 1)
        )
    } else {
        final <- combo
    }

    return(final)
}


#' Plot HoverNet H5AD Segmentation Overlay with Thumbnail
#'
#' @description Creates a side-by-side visualization of the original tissue
#'   thumbnail image and the HoverNet cell segmentation with colored cell type
#'   labels. This function works with HoverNet data stored in h5ad format and
#'   imported as a SpatialExperiment object.
#'
#' @param hovernet A `SpatialExperiment` or `SpatialFeatureExperiment` object
#'   containing HoverNet segmentation data from h5ad format.
#' @param thumbnail_path Path or URL to the thumbnail PNG image file.
#' @param title Optional. Title for the combined plot. If `NULL`, uses
#'   "HoverNet Segmentation". Default is `NULL`.
#' @param point_size Numeric value for the size of points in the segmentation
#'   plot. Default is `0.01`.
#' @param legend_point_size Numeric value for the size of points in the legend.
#'   Default is `2`.
#' @param color_palette Optional. A named vector of colors for cell types. If
#'   `NULL`, uses RColorBrewer palette. Default is `NULL`.
#' @param ncol Number of columns for the side-by-side plots. Default is `2`.
#' @param rel_widths Relative widths of the image and segmentation plots.
#'   Default is `c(1, 1.08)`.
#' @param title_size Font size for the main title. Default is `25`.
#' @param flip_image Logical. Whether to flip the image vertically to match
#'   coordinate system. Default is `TRUE`.
#'
#' @return A `ggplot` object combining the thumbnail image and segmentation
#'   overlay with an optional title.
#'
#' @details The function performs the following steps:
#'   1. Loads the thumbnail PNG image (from local path or URL with caching)
#'   2. Creates a segmentation plot colored by cell type using the `type`
#'      column from `colData`
#'   3. Combines the thumbnail and segmentation side-by-side
#'   4. Adds an optional title
#'
#'   The function expects the `SpatialExperiment` to have:
#'   - `spatialCoords` with columns `x_centroid` and `y_centroid`
#'   - `colData` with a `type` column for cell type classification
#'
#'   Remote PNG files are automatically cached using `BiocFileCache` for
#'   efficient handling without manual downloads.
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual guides
#'   guide_legend theme_void theme coord_fixed element_blank margin
#' @importFrom cowplot plot_grid ggdraw draw_image draw_label get_legend
#' @importFrom SpatialExperiment spatialCoords
#' @importFrom SummarizedExperiment colData
#' @importFrom methods is
#'
#' @examples
#' \dontrun{
#' # From local file
#' thumbnail_path <- "path/to/thumbnail.png"
#' plotHoverNetH5ADOverlay(hn_spe, thumbnail_path = thumbnail_path)
#'
#' # From URL (with automatic caching)
#' thumb_url <- paste0(
#'   "https://store.cancerdatasci.org/hovernet/thumb/",
#'   "TCGA-VG-A8LO-01A-01-DX1.B39A4D64-82A1-4A04-8AB6-918F3058B83B.png"
#' )
#' plotHoverNetH5ADOverlay(hn_spe, thumbnail_path = thumb_url)
#' }
#'
#' @export
plotHoverNetH5ADOverlay <- function(
        hovernet,
        thumbnail_path,
        title = NULL,
        point_size = 0.01,
        legend_point_size = 2,
        color_palette = NULL,
        ncol = 2,
        rel_widths = c(1, 1.08),
        title_size = 25,
        flip_image = TRUE
) {

    # Validate input
    if (!is(hovernet, "SpatialExperiment") &&
        !is(hovernet, "SpatialFeatureExperiment")) {
        stop("'hovernet' must be a SpatialExperiment or ",
             "SpatialFeatureExperiment object.")
    }

    # Check for required columns
    if (!"type" %in% colnames(colData(hovernet))) {
        stop("'hovernet' must contain a 'type' column in colData.")
    }

    coords <- spatialCoords(hovernet)
    if (!all(c("x_centroid", "y_centroid") %in% colnames(coords))) {
        stop("'hovernet' must contain 'x_centroid' and 'y_centroid' in ",
             "spatialCoords.")
    }

    # Check thumbnail path
    if (missing(thumbnail_path) || is.null(thumbnail_path)) {
        stop("'thumbnail_path' is required.")
    }

    # Handle URL or local file for thumbnail
    is_url <- .is_url(thumbnail_path)
    if (is_url)
        thumbnail_path <- .cache_url_file(thumbnail_path)

    if (!file.exists(thumbnail_path))
        stop("Thumbnail file not found: ", thumbnail_path)

    # Prepare data for plotting
    gg <- data.frame(
        x = coords[, "x_centroid"],
        y = coords[, "y_centroid"],
        type = colData(hovernet)$type
    )

    # Determine color palette
    if (is.null(color_palette)) {
        # Get unique types
        unique_types <- unique(gg$type)
        n_types <- length(unique_types)

        # Use RColorBrewer palette
        if (n_types <= 12) {
            color_palette <- RColorBrewer::brewer.pal(
                max(3, n_types), "Paired"
            )[seq_len(n_types)]
        } else {
            # For more than 12 types, cycle through multiple palettes
            color_palette <- grDevices::rainbow(n_types)
        }
        names(color_palette) <- unique_types
    }

    # Create segmentation plot
    p_feat <- ggplot2::ggplot(gg, ggplot2::aes(x, y, col = type)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_manual(
            values = color_palette,
            name = "Cell Type"
        ) +
        ggplot2::guides(
            col = ggplot2::guide_legend(
                override.aes = list(size = legend_point_size)
            )
        ) +
        ggplot2::theme_void()

    # Create version without legend for side-by-side comparison
    p_feat_nolegend <- p_feat +
        ggplot2::theme(legend.position = "none") +
        ggplot2::coord_fixed(ratio = 1)

    # Extract legend
    legend <- cowplot::get_legend(
        p_feat +
            ggplot2::theme(
                legend.position = "right",
                legend.box.margin = ggplot2::margin(0, 0, 0, 10)
            )
    )

    # Prepare image
    img <- magick::image_read(thumbnail_path)
    if (flip_image) {
        img <- magick::image_flip(img)
    }

    # Create image plot
    p_img <- cowplot::ggdraw() + cowplot::draw_image(img)

    # Combine plots side by side
    combo_plots <- cowplot::plot_grid(
        p_img, p_feat_nolegend,
        ncol = ncol,
        align = "hv",
        axis = "tblr",
        rel_widths = rel_widths
    )

    # Add legend on the right
    combo <- cowplot::plot_grid(
        combo_plots, legend,
        ncol = 2,
        rel_widths = c(1, 0.15)
    )

    # Add title if provided
    if (is.null(title)) {
        title <- "HoverNet Segmentation"
    }

    if (!is.null(title) && nchar(title) > 0) {
        title_plot <- cowplot::ggdraw() +
            cowplot::draw_label(
                title,
                fontface = "bold",
                x = 0.5,
                hjust = 0.5,
                size = title_size
            )

        final <- cowplot::plot_grid(
            title_plot, combo,
            ncol = 1,
            rel_heights = c(0.08, 1)
        )
    } else {
        final <- combo
    }

    return(final)
}
