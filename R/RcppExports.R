# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

assign2visium_spots <- function(mole_coords, array_coords, img_coords, spot_radius, thread_num = 1L, verbose = TRUE) {
    .Call(`_blurrr_assign2visium_spots`, mole_coords, array_coords, img_coords, spot_radius, thread_num, verbose)
}

assign2visium_subspots <- function(mole_coords, mole_assigned_barcodes, subspot_ids, img_subspot_coords, subspot_assigned_barcodes, img_spot_coords, spot_barcodes, spot_radius, thread_num, verbose) {
    .Call(`_blurrr_assign2visium_subspots`, mole_coords, mole_assigned_barcodes, subspot_ids, img_subspot_coords, subspot_assigned_barcodes, img_spot_coords, spot_barcodes, spot_radius, thread_num, verbose)
}

bin_transcripts <- function(unique_spot_names, unique_gene_names, assignment_mole_names, assignment_spot_names, feature_mole_names, feature_gene_names, thread_num = 1L, verbose = TRUE) {
    .Call(`_blurrr_bin_transcripts`, unique_spot_names, unique_gene_names, assignment_mole_names, assignment_spot_names, feature_mole_names, feature_gene_names, thread_num, verbose)
}

