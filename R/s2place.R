#' Convert place name to location on sphere
#'
#' @param placename What to find.
#' @param online Logical. Should an online lookup be performed if needed?
#' @param verbose Logical. Should verbose output be produced?
#'
#' @export
s2place <- function(placename, online = TRUE, verbose = FALSE) {
  stopifnot(is.character(placename))
  pn <- tolower(placename)
  rslt <- try(globe::place(pn), silent = TRUE)
  if(!inherits(rslt, "try-error"))
    return(rslt)
  if(online && requireNamespace("RgoogleMaps", quietly = TRUE)){
    if(verbose){
      message("Didn't find ", sQuote(placename),
              " among builtin placenames in globe package -- trying online via Google Maps API...")
    }
    rslt <- RgoogleMaps::getGeoCode(pn, verbose = verbose)
    if(!anyNA(rslt)) return(list(lon = rslt["lon"], lat = rslt["lat"]))
  } else{
    if(verbose){
      message("Install package RgoogleMaps to look for ", sQuote(placename), " online")
    }
  }
  stop(paste("Unrecognised placename", sQuote(placename),
             "-- check globe package for available builtin places"),
       call.=FALSE)
}
