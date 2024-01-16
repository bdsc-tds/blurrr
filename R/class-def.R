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
  "XeniumMolecule",
  slots = c(
    pos = "data.frame",
    rotMtx = "NullOrMatrix",
    alignedFile = "NullOrCharacter"
  )
)
setClassUnion("NullOrXeniumMolecule", c("NULL", "XeniumMolecule"))


#' @export
#' 
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
setClass(
  "XeniumCell",
  slots = c(
    sce = "SingleCellExperiment"
  ),
  contains = c("XeniumMolecule")
)
setClassUnion("NullOrXeniumCell", c("NULL", "XeniumCell"))


#' @export
setClass(
  "AssignedXeniumMolecule",
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
setClassUnion("NullOrAssignedXeniumMolecule", c("NULL", "AssignedXeniumMolecule"))


#' @export
setClass(
  "BinXenium",
  slots = c(
    vsInfo = "VisiumInfo",
    xnCell = "NullOrXeniumCell",
    xnTrans = "NullOrXeniumMolecule",
    is2Subspot = "logical",
    assignedXnCell = "NullOrAssignedXeniumMolecule",
    assignedXnTrans = "NullOrAssignedXeniumMolecule"
  )
)
