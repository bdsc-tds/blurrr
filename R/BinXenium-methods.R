# getters ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_vsInfo",
  "BinXenium",
  function(x) x@vsInfo
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnCell",
  "BinXenium",
  function(x) x@xnCell
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnTrans",
  "BinXenium",
  function(x) x@xnTrans
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_is2Subspot",
  "BinXenium",
  function(x) x@is2Subspot
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignedXnCell",
  "BinXenium",
  function(x) x@assignedXnCell
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignedXnTrans",
  "BinXenium",
  function(x) x@assignedXnTrans
)


# others ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "assign2visium",
  c("BinXenium", "numeric", "logical", "logical"),
  function(x, cores, force, verbose) {
    if (!is.null(get_xnCell(x))) {
      if (verbose) {
        message("Assigning cells...")
      }
      
      x@assignedXnCell <- .assign2visium(
        assigned.xn = get_assignedXnCell(x),
        vs.pos = get_vsPos(x),
        xn.pos = get_xnPos(x, is.cell = TRUE),
        is2subspot = get_is2Subspot(x),
        subspot.pos = get_subspotPos(x),
        is.cell = TRUE,
        spot_radius = get_scalef(x, "d") / 2,
        cores = cores,
        force = force,
        verbose = verbose
      )
    }
    
    if (!is.null(get_xnTrans(x))) {
      if (verbose) {
        message("Assigning transcripts...")
      }
      
      x@assignedXnTrans <- .assign2visium(
        assigned.xn = get_assignedXnTrans(x),
        vs.pos = get_vsPos(x),
        xn.pos = get_xnPos(x, is.cell = FALSE),
        is2subspot = get_is2Subspot(x),
        subspot.pos = get_subspotPos(x),
        is.cell = FALSE,
        spot_radius = get_scalef(x, "d") / 2,
        cores = cores,
        force = force,
        verbose = verbose
      )
    }
    
    x
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "assign2visium",
  c("BinXenium", "numeric", "missing", "missing"),
  function(x, cores, force = FALSE, verbose = TRUE) {
    assign2visium(
      x = x,
      cores = cores,
      force = FALSE,
      verbose = TRUE
    )
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "assign2visium",
  c("BinXenium", "numeric", "logical", "missing"),
  function(x, cores, force, verbose = TRUE) {
    assign2visium(
      x = x,
      cores = cores,
      force = force,
      verbose = TRUE
    )
  }
)

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "assign2visium",
  c("BinXenium", "numeric", "missing", "logical"),
  function(x, cores, force = FALSE, verbose) {
    assign2visium(
      x = x,
      cores = cores,
      force = FALSE,
      verbose = verbose
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "numeric", "logical", "logical"),
  function(x, dirname, cores, cache, verbose) {
    dir.create(file.path(dirname), showWarnings = FALSE)
    
    if (cache) {
      saveRDS(
        object = x,
        file = file.path(dirname, "bin_xenium.rds")
      )
    }
    
    if (!is.null(get_xnCell(x)) && !is.null(get_assignedXnCell(x))) {
      .save2disk_cell(
        xn.sce = get_xn_sce(x),
        assigned.xn.cell = get_assignedXnCell(x),
        dirname = dirname,
        is.subspot = get_is2Subspot(x),
        cores = cores,
        verbose = verbose
      )
    }
    
    if (!is.null(get_xnTrans(x)) && !is.null(get_assignedXnTrans(x))) {
      NULL
    }
    
    invisible(NULL)
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "numeric", "logical", "missing"),
  function(x, dirname, cores, cache, verbose = TRUE) {
    save2disk(
      x = x,
      dirname = dirname,
      cores = cores,
      cache = cache,
      verbose = TRUE
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "numeric", "missing", "missing"),
  function(x, dirname, cores, cache = FALSE, verbose = TRUE) {
    save2disk(
      x = x,
      dirname = dirname,
      cores = cores,
      cache = FALSE,
      verbose = TRUE
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "missing", "logical", "missing"),
  function(x, dirname, cores = 1, cache, verbose = TRUE) {
    save2disk(
      x = x,
      dirname = dirname,
      cores = 1,
      cache = cache,
      verbose = TRUE
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "missing", "missing", "missing"),
  function(x, dirname, cores = 1, cache = FALSE, verbose = TRUE) {
    save2disk(
      x = x,
      dirname = dirname,
      cores = 1,
      cache = FALSE,
      verbose = TRUE
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "missing", "logical", "numeric"),
  function(x, dirname, cores = 1, cache, verbose) {
    save2disk(
      x = x,
      dirname = dirname,
      cores = 1,
      cache = cache,
      verbose = verbose
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "save2disk",
  c("BinXenium", "character", "missing", "missing", "numeric"),
  function(x, dirname, cores = 1, cache = FALSE, verbose) {
    save2disk(
      x = x,
      dirname = dirname,
      cores = 1,
      cache = cache,
      verbose = verbose
    )
  }
)
