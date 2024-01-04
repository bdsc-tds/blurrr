setClass(
  "VisiumInfo",
  slots = c(
    pos = "data.frame",
    scalef = "list"
  )
)


setClass(
  "XeniumMolecular",
  slots = c(
    pos = "data.frame"
  )
)
setClassUnion("NullOrXeniumMolecular", c("NULL", "XeniumMolecular"))


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
    toSpot = "logical",
    toSubspot = "logical"
  )
)
