#' @name BinXenium-methods
#' 
#' @title Methods for class `BinXenium`
#' 
#' @aliases
#' get_vsInfo
#' get_xnCell
#' get_xnTrans
#' get_is2Subspot
#' get_assignedXnCell
#' get_assignedXnTrans
#' assign2visium
#' subset2assigned
#' plot_mole
#' save2disk
#' 
#' @description 
#' 
NULL

# getters ---------------------------------------------

#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_vsInfo",
  "BinXenium",
  function(x) x@vsInfo
)


#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnCell",
  "BinXenium",
  function(x) x@xnCell
)


#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_xnTrans",
  "BinXenium",
  function(x) x@xnTrans
)


#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_is2Subspot",
  "BinXenium",
  function(x) x@is2Subspot
)


#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignedXnCell",
  "BinXenium",
  function(x) x@assignedXnCell
)


#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_assignedXnTrans",
  "BinXenium",
  function(x) x@assignedXnTrans
)


# assign2visium ---------------------------------------------

#' @rdname BinXenium-methods
#'
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


# subset2assigned ---------------------------------------------

#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom magrittr `%>%`
#' @importFrom dplyr right_join select
setMethod(
  "subset2assigned",
  c("BinXenium", "character", "character"),
  function(
    x,
    mole = c("cell", "trans"),
    level = c("spot", "subspot")
  ) {
    mole = match.arg(mole)
    level = match.arg(level)
    
    if (mole == "cell") {
      assert_that(!is.null(get_xnCell(x)))
      
      if (level == "spot") {
        assert_that(
          !is.null(get_assignment2Spots(x, TRUE))
        )
        
        .filtered <- get_assignment2Spots(x, TRUE)
      } else {
        assert_that(
          get_is2Subspot(x),
          !is.null(get_assignment2Subspots(x, TRUE))
        )
        
        .filtered <- get_assignment2Subspots(x, TRUE)
      }
      
      is.cell <- TRUE
      .id <- "cell_id"
    } else {
      assert_that(!is.null(get_xnTrans(x)))
      
      if (level == "spot") {
        assert_that(
          !is.null(get_assignment2Spots(x, FALSE))
        )
        
        .filtered <- get_assignment2Spots(x, FALSE)
      } else {
        assert_that(
          get_is2Subspot(x),
          !is.null(get_assignment2Subspots(x, FALSE))
        )
        
        .filtered <- get_assignment2Subspots(x, FALSE)
      }
      
      is.cell <- FALSE
      .id <- "transcript_id"
    }
    
    new(
      "XeniumMolecule",
      pos = get_xnPos(x, is.cell) %>%
        right_join(
          .filtered %>%
            select(id),
          by = setNames(nm = .id, "id")
        ),
      rotMtx = NULL,
      alignedFile = NULL
    )
  }
)


# plot_mole ---------------------------------------------

#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
#'
#' @importFrom assertthat assert_that
#' @importFrom ggpubr ggarrange
setMethod(
  "plot_mole",
  c("BinXenium", "character", "character", "character"),
  function(
    x,
    mole = c("cell", "trans"),
    mode = c("raw", "assigned", "both"),
    res = c("fullres", "hires", "lowres")
  ) {
    mole = match.arg(mole)
    mode = match.arg(mode)
    res = match.arg(res)
    
    if (mole == "cell" && mode %in% c("raw", "both")) {
      assert_that(!is.null(get_xnCell(x)))
      is.cell <- TRUE
      .raw.data <- get_xnCell(x)
    }
    
    if (mole == "cell" && mode %in% c("assigned", "both")) {
      assert_that(!is.null(get_assignedXnCell(x)))
      is.cell <- TRUE
      .assigned.data <- get_assignedXnCell(x)
    }
    
    if (mole == "trans" && mode %in% c("raw", "both")) {
      assert_that(!is.null(get_xnTrans(x)))
      is.cell <- FALSE
      .raw.data <- get_xnTrans(x)
    }
    
    if (mole == "trans" && mode %in% c("assigned", "both")) {
      assert_that(!is.null(get_assignedXnTrans(x)))
      is.cell <- FALSE
      .assigned.data <- get_assignedXnTrans(x)
    }
    
    # plot background
    .bg <- .plot_bg(
      vs.info = get_vsInfo(x),
      res = res
    )
    
    raw.plot <- NULL
    if (mode %in% c("raw", "both")) {
      raw.plot <- .bg + .plot_mole(
        mole.data = get_xnPos(x, is.cell),
        is.cell = is.cell,
        res = res
      )
      
      if (mode == "raw") {
        return(raw.plot)
      }
    }
    
    assigned.plot <- NULL
    if (mode %in% c("assigned", "both")) {
      assigned.plot <- .bg + .plot_mole(
        mole.data = get_xnPos(
          subset2assigned(
            x = x,
            mole = mole,
            level = "spot"
          )
        ),
        is.cell = is.cell,
        res = res
      )
      
      if (mode == "raw") {
        return(assigned.plot)
      }
    }
    
    return(
      ggarrange(
        raw.plot,
        assigned.plot,
        ncol = 2
      )
    )
  }
)


# save2disk ---------------------------------------------

#' @rdname BinXenium-methods
#'
#' @include generic-def.R class-def.R
#' 
#' @export
#' 
#' @importFrom magrittr `%>%`
#' @importFrom dplyr select
setMethod(
  "save2disk",
  c("BinXenium", "character", "numeric", "logical", "logical"),
  function(x, dirname, cores, cache, verbose) {
    dir.create(file.path(dirname), showWarnings = FALSE)
    dirname <- normalizePath(dirname)
    
    if (cache) {
      saveRDS(
        object = x,
        file = file.path(dirname, "bin_xenium.rds")
      )
    }
    
    if (!is.null(get_xnCell(x)) && !is.null(get_assignedXnCell(x))) {
      .save2disk(
        assigned.xn = get_assignedXnCell(x),
        dirname = dirname,
        is.cell = TRUE,
        is.subspot = get_is2Subspot(x),
        xn.sce = get_xn_sce(x),
        cores = cores,
        verbose = verbose
      )
    }
    
    if (!is.null(get_xnTrans(x)) && !is.null(get_assignedXnTrans(x))) {
      xn.pos = get_xnPos(x, is.cell = FALSE) %>%
        select(
          transcript_id,
          feature_name
        )
      
      .save2disk(
        assigned.xn = get_assignedXnTrans(x),
        dirname = dirname,
        is.cell = FALSE,
        is.subspot = get_is2Subspot(x),
        xn.pos = xn.pos,
        cores = cores,
        verbose = verbose
      )
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
      cache = FALSE,
      verbose = verbose
    )
  }
)
