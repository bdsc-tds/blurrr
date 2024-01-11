setClassUnion("NullOrMatrix", c("NULL", "matrix"))
setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrDataFrame", c("NULL", "data.frame"))


#' @export
setClass(
  "VisiumInfo",
  slots = c(
    pos = "data.frame",
    arrayDirec = "list",
    imgRad = "numeric",
    scalef = "list",
    subspotPos = "NullOrDataFrame"
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
  "AssignedXeniumMolecular",
  slots = c(
    assignment2Spots = "NullOrDataFrame",
    ambiAssignment2Spots = "NullOrDataFrame",
    countOfSpots = "NullOrDataFrame",
    assignment2Subspots = "NullOrDataFrame",
    ambiAssignment2Subspots = "NullOrDataFrame",
    countOfSubspots = "NullOrDataFrame",
    propAssigned = "numeric"
  )
)
setClassUnion("NullOrAssignedXeniumMolecular", c("NULL", "AssignedXeniumMolecular"))


#' @export
setClass(
  "BinXenium",
  slots = c(
    vsInfo = "VisiumInfo",
    xnCell = "NullOrXeniumCell",
    xnTrans = "NullOrXeniumMolecular",
    is2Subspot = "logical",
    assignedXnCell = "NullOrAssignedXeniumMolecular",
    assignedXnTrans = "NullOrAssignedXeniumMolecular"
  )
)