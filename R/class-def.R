setClassUnion("NullOrMatrix", c("NULL", "matrix"))
setClassUnion("NullOrCharacter", c("NULL", "character"))
setClassUnion("NullOrDataFrame", c("NULL", "data.frame"))


#' @export
setClass(
  "Template",
  slots = c(
    type = "character",
    pos = "data.frame",
    arrayDirec = "list",
    imgRad = "numeric",
    scalef = "list",
    subspotPos = "NullOrDataFrame"
  )
)


#' @export
setClass(
  "Molecule",
  slots = c(
    type = "character",
    sce = "SingleCellExperiment",
    rotMtx = "NullOrMatrix",
    alignedFile = "NullOrCharacter"
  )
)
setClassUnion("NullOrMolecule", c("NULL", "Molecule"))


#' @export
setClass(
  "MappedMolecule",
  slots = c(
    map = "NullOrDataFrame",
    mapInDoubt = "NullOrDataFrame",
    cnt = "NullOrDataFrame",
    propAssigned = "numeric"
  )
)
setClassUnion("NullOrMappedMolecule", c("NULL", "MappedMolecule"))


#' @export
setClass(
  "Blurrr",
  slots = c(
    template = "Template",
    molecule = "Molecule",
    is2Subspot = "logical",
    mappedMolecule = "NullOrMappedMolecule",
    mappedMolecule2subspot = "NullOrMappedMolecule"
  )
)
