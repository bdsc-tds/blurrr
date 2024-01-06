# For `VisiumInfo` ---------------------------------------------

#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_arrayDirec",
  c("VisiumInfo", "character"),
  function(x, key = c("row", "col", "both")) {
    key = match.arg(key)
    
    switch(
      key,
      "row" = {
        return(x@arrayDirec$row)
      },
      "col" = {
        return(x@arrayDirec$col)
      },
      {
        return(x@arrayDirec)
      }
    )
  }
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_imgRad",
  "VisiumInfo",
  function(x) x@imgRad
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_scalef",
  c("VisiumInfo", "character"),
  function(x, key = c("hires", "lowres", "diameter")) {
    key = match.arg(key)
    
    switch(
      key,
      "hires" = {
        return(x@scalef$tissue_hires_scalef)
      },
      "lowres" = {
        return(x@scalef$tissue_lowres_scalef)
      },
      "diameter" = {
        return(x@scalef$spot_diameter_fullres)
      },
      {
        return(NULL)
      }
    )
  }
)


# For `BinXenium` ---------------------------------------------


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_arrayDirec",
  c("BinXenium", "character"),
  function(x, key) get_arrayDirec(vsInfo(x), key)
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_imgRad",
  "BinXenium",
  function(x) get_imgRad(vsInfo(x))
)


#' @include generic-def.R class-def.R
#' 
#' @export
setMethod(
  "get_scalef",
  c("BinXenium", "character"),
  function(x, key) get_scalef(vsInfo(x), key)
)
