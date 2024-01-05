setClassUnion("NullOrMatrix", c("NULL", "matrix"))
setClassUnion("NullOrCharacter", c("NULL", "character"))


#' @export
setClass(
  "VisiumInfo",
  slots = c(
    pos = "data.frame",
    arrayDirec = "list",
    imgRad = "numeric",
    scalef = "list"
  )
)


#' @export
setClass(
  "XeniumMolecular",
  slots = c(
    pos = "data.frame",
    rotMtx = "NullOrMatrix",
    alignedFile = "NullOrCharacter"
  )
)
setClassUnion("NullOrXeniumMolecular", c("NULL", "XeniumMolecular"))


#' @export
#' 
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setClass(
  "XeniumCell",
  slots = c(
    sce = "SingleCellExperiment"
  ),
  contains = c("XeniumMolecular")
)
setClassUnion("NullOrXeniumCell", c("NULL", "XeniumCell"))


#' @export
setClass(
  "BinXenium",
  slots = c(
    vsInfo = "VisiumInfo",
    xnCell = "NullOrXeniumCell",
    xnTrans = "NullOrXeniumMolecular",
    is2Spot = "logical",
    is2Subspot = "logical"
  )
)
