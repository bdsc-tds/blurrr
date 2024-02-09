#' @keywords internal
#' 
#' @importFrom magrittr `%>%`
#' @importFrom dplyr rename mutate left_join select
.assign2template <- function(
    assigned.xn,
    vs.pos,
    xn.pos,
    is2subspot,
    subspot.pos,
    is.cell,
    spot_radius,
    cores,
    force = FALSE,
    verbose = TRUE
) {
  mole.id <- ifelse(
    is.cell,
    "cell_id",
    "transcript_id"
  )
  
  to.spot <- FALSE
  if (force || is.null(assigned.xn) || is.null(get_assignment2Spots(assigned.xn)) || is.null(get_countOfSpots(assigned.xn))) {
    to.spot <- TRUE
  }
  
  to.subspot <- FALSE
  if (is2subspot && (force || to.spot || is.null(get_assignment2Subspots(assigned.xn)) || is.null(get_countOfSubspots(assigned.xn)))) {
    to.subspot <- TRUE
  }
  
  binned.subspot.avail <- TRUE
  if (is.null(assigned.xn)) {
    binned.subspot.avail <- FALSE
  }
  
  if (to.spot) {
    spot_ids = vs.pos$barcode
    
    .assigned2spots <- assign2visium_spots(
      mole_coords = as.matrix(
        xn.pos[get_coords_names(
          is.xn = TRUE,
          is.cell = is.cell,
          prefix = "array_aligned",
          use.names = "f"
        )[[1]]]
      ),
      array_coords = as.matrix(
        vs.pos[c("array_col", "array_row")]
      ),
      img_coords = as.matrix(
        vs.pos[get_coords_names(
          is.xn = FALSE,
          prefix = "array_aligned",
          use.names = "f"
        )[[1]]]
      ),
      spot_radius = spot_radius,
      thread_num = cores,
      verbose = verbose
    )
    
    .assignment2Spots <- NULL
    if (dim(.assigned2spots$assignment2Spots)[1] > 0) {
      .assignment2Spots <- as.data.frame(.assigned2spots$assignment2Spots) %>%
        rename(
          id = V1,
          barcode = V2
        ) %>%
        mutate(
          id = xn.pos[[mole.id]][id],
          barcode = spot_ids[barcode]
        )
    }
    
    .ambiAssignment2Spots <- NULL
    if (dim(.assigned2spots$ambiAssignment2Spots)[1] > 0) {
      .ambiAssignment2Spots <- as.data.frame(.assigned2spots$ambiAssignment2Spots) %>%
        rename(
          id = V1,
          barcode = V2
        ) %>%
        mutate(
          id = xn.pos[["mole.id"]][id],
          barcode = spot_ids[barcode]
        )
    }
    
    .countOfSpots <- NULL
    if (dim(.assigned2spots$countOfSpots)[1] > 0) {
      .countOfSpots <- as.data.frame(.assigned2spots$countOfSpots) %>%
        rename(
          barcode = V1,
          count = V2
        ) %>%
        mutate(
          barcode = spot_ids[barcode]
        )
    }
    
    assigned2spots <- list(
      assignment2Spots = .assignment2Spots,
      ambiAssignment2Spots = .ambiAssignment2Spots,
      countOfSpots = .countOfSpots,
      propAssigned = .assigned2spots$propAssigned
    )
  } else {
    assigned2spots <- list(
      assignment2Spots = get_assignment2Spots(assigned.xn),
      ambiAssignment2Spots = get_ambiAssignment2Spots(assigned.xn),
      countOfSpots = get_countOfSpots(assigned.xn),
      propAssigned = get_propAssigned(assigned.xn)
    )
  }
  
  if (to.subspot) {
    .assigned2subspots <- assign2visium_subspots(
      mole_coords = as.matrix(
        left_join(
          assigned2spots$assignment2Spots,
          xn.pos[c(
            mole.id,
            get_coords_names(
              is.xn = TRUE,
              is.cell = is.cell,
              use.names = "f"
            )[[1]]
          )],
          by = c("id" = mole.id)
        )[get_coords_names(
          is.xn = TRUE,
          is.cell = is.cell,
          use.names = "f"
        )[[1]]]
      ),
      mole_assigned_barcodes = assigned2spots$assignment2Spots[["barcode"]],
      subspot_ids = as.matrix(
        subspot.pos[c("subspot_id", "third_pt_id")]
      ),
      img_subspot_coords = as.matrix(
        subspot.pos[get_coords_names(
          is.xn = FALSE,
          is.subspot = TRUE
        )[[1]]]
      ),
      subspot_assigned_barcodes = subspot.pos[["barcode"]],
      img_spot_coords = as.matrix(
        vs.pos[get_coords_names(
          is.xn = FALSE,
          use.names = "f"
        )[[1]]]
      ),
      spot_barcodes = vs.pos[["barcode"]],
      spot_radius = spot_radius,
      thread_num = cores,
      verbose = verbose
    )
    
    .assignment2Subspots <- NULL
    if (dim(.assigned2subspots$assignment2Subspots)[1] > 0) {
      .assignment2Subspots <- as.data.frame(.assigned2subspots$assignment2Subspots) %>%
        rename(
          id = V1,
          subspot_barcode_idx = V2
        ) %>%
        mutate(
          id = assigned2spots$assignment2Spots$id[id],
          subspot_barcode = subspot.pos$subspot_barcode[subspot_barcode_idx],
          subspot_barcode_bayesspace = subspot.pos$subspot_barcode_bayesspace[subspot_barcode_idx]
        ) %>%
        select(
          -subspot_barcode_idx
        )
    }
    
    .ambiAssignment2Subspots <- NULL
    if (dim(.assigned2subspots$ambiAssignment2Subspots)[1] > 0) {
      .ambiAssignment2Subspots <- as.data.frame(.assigned2subspots$ambiAssignment2Subspots) %>%
        rename(
          id = V1,
          subspot_barcode_idx = V2
        ) %>%
        mutate(
          id = assigned2spots$assignment2Spots$id[id],
          subspot_barcode = subspot.pos$subspot_barcode[subspot_barcode_idx],
          subspot_barcode_bayesspace = subspot.pos$subspot_barcode_bayesspace[subspot_barcode_idx]
        ) %>%
        select(
          -subspot_barcode_idx
        )
    }
    
    .countOfSubspots <- NULL
    if (dim(.assigned2subspots$countOfSubspots)[1] > 0) {
      .countOfSubspots <- as.data.frame(.assigned2subspots$countOfSubspots) %>%
        rename(
          subspot_barcode_idx = V1,
          count = V2
        ) %>%
        mutate(
          subspot_barcode = subspot.pos$subspot_barcode[subspot_barcode_idx],
          subspot_barcode_bayesspace = subspot.pos$subspot_barcode_bayesspace[subspot_barcode_idx]
        ) %>%
        select(
          -subspot_barcode_idx
        ) %>%
        select(
          subspot_barcode,
          subspot_barcode_bayesspace,
          count
        )
    }
    
    assigned2subspots <- list(
      assignment2Subspots = .assignment2Subspots,
      ambiAssignment2Subspots = .ambiAssignment2Subspots,
      countOfSubspots = .countOfSubspots
    )
  } else if (binned.subspot.avail) {
    assigned2subspots <- list(
      assignment2Subspots = get_assignment2Subspots(assigned.xn),
      ambiAssignment2Subspots = get_ambiAssignment2Subspots(assigned.xn),
      countOfSubspots = get_countOfSubspots(assigned.xn)
    )
  } else {
    assigned2subspots <- list(
      assignment2Subspots = NULL,
      ambiAssignment2Subspots = NULL,
      countOfSubspots = NULL
    )
  }
  
  new(
    "AssignedXeniumMolecule",
    assignment2Spots = assigned2spots$assignment2Spots,
    ambiAssignment2Spots = assigned2spots$ambiAssignment2Spots,
    countOfSpots = assigned2spots$countOfSpots,
    assignment2Subspots = assigned2subspots$assignment2Subspots,
    ambiAssignment2Subspots = assigned2subspots$ambiAssignment2Subspots,
    countOfSubspots = assigned2subspots$countOfSubspots,
    propAssigned = assigned2spots$propAssigned
  )
}
