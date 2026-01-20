#' Import CONCH as a SpatialFeatureExperiment
#'
#' @description This function reads CONCH output stored in an HDF5 (`.h5`) file
#'   and converts it into a `SpatialFeatureExperiment` object for downstream
#'   analysis and visualization. The function extracts spatial coordinates and
#'   feature embeddings from the file, ensuring proper formatting for the SFE
#'   class.
#'
#' @param file_path `character(1)` Path to the .h5 file containing CONCH
#'   output.
#' @param patch_size `numeric(1)` The width/height of the patch in pixels
#'   (default 224). Used to set the spatial metadata.
#'
#' @importFrom S4Vectors metadata<-
#'
#' @return A `SpatialFeatureExperiment` class object.
#'
#' @export
CONCH <- function(file_path, patch_size = 224) {
    checkInstalled("rhdf5")
    # 1. Read Data
    # h5read will read these based on the dimensions in the file
    # coords: usually 2 x N (x, y)
    # features: usually 768 x N (embedding_dim x patches)
    coords_raw <- rhdf5::h5read(file_path, "coords")
    
    # 2. Prepare Spatial Coordinates
    # SFE expects spatialCoords to be (N_samples x 2)
    # If input is (2 x N), we transpose it.
    scFUN <- if (identical(nrow(coords_raw), 2L)) t else identity
    spatial_coords <- scFUN(coords_raw)
    colnames(spatial_coords) <- c("x", "y")
    
    # 3. Prepare Assay Data (Features)
    # SFE assays expect (N_features x N_samples)
    # Your input is already 768 x 3692, so it likely needs no transposition.
    assay_data <- rhdf5::h5read(file_path, "features")
    
    # Safety check: if dimensions were flipped during generation
    if (!identical(ncol(assay_data), nrow(spatial_coords))) {
        warning(
            "Dimensions mismatch detected. ",
            "Attempting to transpose features matrix..."
        )
        assay_data <- t(assay_data)
    }
    
    # Assign generic names to avoid empty rownames issues later
    rownames(assay_data) <- paste0("dim_", seq_len(nrow(assay_data)))
    colnames(assay_data) <- paste0("patch_", seq_len(ncol(assay_data)))
    
    checkInstalled("SpatialFeatureExperiment")
    # 4. Construct the SpatialFeatureExperiment Object
    sfe <- SpatialFeatureExperiment::SpatialFeatureExperiment(
        assays = list(conch = assay_data),
        spatialCoords = spatial_coords,
        # Optional: Set unit if known, usually 'px' for raw slides
        unit = "px" 
    )
    
    # 5. Add Metadata (Optional but recommended for reproducibility)
    metadata(sfe)$source_file <- file_path
    metadata(sfe)$patch_size <- patch_size
    metadata(sfe)$model <- "CONCH"
    
    sfe
}