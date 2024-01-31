#' @name Blurrr
#' 
#' @title Create a `Blurrr` object.
#'
#' @param mole The type of molecules to blur. Possible values include
#'   `xenium_cell`, `xenium_transcript` and `visium_hd`, which are mapped via
#'   `match.arg`.
#' @param blur The type of blurring. Possible values include `visium` and
#'   `visium_hd`, which are mapped via `match.arg`.
#' @param mole.dir The directory to the data of molecules. For any \code{mole}
#'   it should be the deepest directory containing all the generated data.
#' @param blur.template.dir The directory to the template blurring. Any provided
#'   path indicates that the molecules to blur have been aligned to this
#'   template, and it should be the deepest directory containing all the
#'   generated data. Otherwise a provided `NULL` indicates that the molecules
#'   will be blurred following specific instructions (see details).
#' @param blur2subspot Whether to blur to subspot level or not. If yes and
#'   \code{blur} is `visium_hd`, set other arguments to control the type of
#'   subspots (see details).
#' @param ... Extra arguments (see details).
#' 
#' @details
#' Valid extra arguments passed through \code{...} are different dependent on
#' values passed to other arguments.
#' 
#' @returns A Blurrr object.
NULL


#' @export
#' 
#' @importFrom assertthat assert_that
Blurrr <- function(
    mole = c("xenium_cell", "xenium_transcript", "visium_hd"),
    blur = c("visium", "visium_hd"),
    mole.dir = ".",
    blur.template.dir = NULL,
    blur2subspot = FALSE,
    ...
) {
  .args <- list(...)
  mole <- match.arg(mole)
  blur <- match.arg(blur)
  
  assert_that(
    file.exists(mole.dir),
    is.null(blur.template.dir) || file.exists(blur.template.dir),
    is.logical(blur2subspot)
  )
  
  # load or create a blurring template
  if (is.null(blur.template.dir)) {
    .template <- .create_blur_template()
  } else {
    .template <- do.call(
      .load_blur_template,
      c(
        list(
          dirname = blur.template.dir,
          blur = blur,
          blur2subspot = blur2subspot
        )
      )
    )
  }
  
  # load molecules
  .molecule <- .load_molecules()
  
  # constructor
  new(
    "Blurrr",
    template = .template,
    molecule = .molecule,
    is2Subspot = blur2subspot,
    mappedMolecule = NULL,
    mappedMolecule2subspot = NULL
  )
}
