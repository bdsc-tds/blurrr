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
#'   generated data. For instance, for \code{blur = visium},
#'   \code{blur.template.dir} should end with `outs`, while for
#'   \code{blur = visium_hd}, \code{blur.template.dir} should end with
#'   `square*`. Otherwise a provided `NULL` indicates that the molecules will
#'   be blurred following specific instructions (see details).
#' @param blur2subspot Whether to blur to subspot level or not. If yes and
#'   \code{blur = visium_hd}, set other arguments to control the type of
#'   subspots (see details).
#' @param rm.mole.feats.pat A vector of patterns to remove for the features of
#'   molecules. Set to `NULL` if everything is to be kept.
#' @param ... Extra arguments (see details).
#' 
#' @details
#' Valid extra arguments passed through \code{...} are listed below.
#'   
#'   1. General arguments.
#'   \itemize{
#'     \item{\code{dbg.load.mole.nrows}}{
#'       For debugging: the number of rows of the molecule data to load.
#'     }
#'   }
#'   
#'   2. Regarding the alignment of molecules to a blurring template, there are
#'   two ways to specify such alignment, whose arguments are mutually exclusive.
#'   
#'   If the molecules are aligned to the blurring template using a rotation
#'   matrix:
#'   \itemize{
#'     \item{\code{rot.mtx}}{}
#'     \item{\code{rot.mtx.to = c("fullres", "hires", "lowres")}}{
#'       The resolution of the image in the blurring template that the molecules
#'       are aligned to.
#'     }
#'   }
#'   
#'   If the molecules are aligned to the blurring template manually:
#'   \itemize{
#'     \item{\code{aligned.file}}{}
#'     \item{\code{aligned.coords.names = c(x = , y = )}}{}
#'     \item{\code{aligned.to = c("fullres", "hires", "lowres")}}{
#'       The resolution of the image in the blurring template that the molecules
#'       are aligned to.
#'     }
#'   }
#'   
#'   When in-place blurring of molecules is of interest, such as the case where
#'   \code{blur.template.dir} is `NULL`, users can leave the above arguments
#'   undefined.
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
    rm.mole.feats.pat = c("^NegControl.*", "^BLANK.*"),
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
  .molecule <- do.call(
    .load_molecules,
    c(
      list(
        dirname = mole.dir,
        image.rad = get_imgRad(.template),
        scalef = list(
          hires = get_scalef(.template, "h"),
          lowres = get_scalef(.template, "l")
        ),
        mole = mole,
        rm.mole.feats.pat = rm.mole.feats.pat
      ),
      .parse_extra_args_mole(args)
    )
  )
  
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
