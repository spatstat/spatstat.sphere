#' Make s2index for points on the sphere (beta version)
#'
#' Create a s2index for a collection of points on the sphere. This is a beta
#' version of a spatial index that stores the s2cellid at various levels for a
#' given collection of points. Expect changes to the API or even complete removal
#' of the function. Only for testing at the moment!
#'
#' @param x Points on the sphere in a standard format.
#' @param min_level Integer between 0 and 30 (incl). Must be smaller than or
#' equal to `max_level`.
#' @param max_level Integer between 0 and 30 (incl). Must be larger than or
#' equal to `min_level`.
#'
#' @return An environment with additional class `s2index`.
#' @export
s2index <- function(x, min_level = 0L, max_level = 30L){
  if(inherits(x, "s2index"))
    return(x)
  stopifnot(is.numeric(c(min_level, max_level)) &&
            length(c(min_level, max_level)) == 2L &&
            min_level >= 0L && min_level <= max_level & max_level <= 30L)
  x <- make_s2coords(x)
  lev <- seq(min_level, max_level)
  indexlist <- list()
  ## For loop over each level:
  for(l in lev){
    ids <- s2cellid(x, l)$id
    ## Indices where each unique cellid occurs:
    indices <- split(seq_along(ids), ids)
    indexlist <- c(indexlist, indices)
  }
  env <- list2env(indexlist)
  class(env) <- c("s2index", class(env))
  return(env)
}

#' Print object of class s2index
#'
#' Print object of class s2index.
#'
#' @param x Object of class `s2index`.
#' @param ... Ignored.
#'
#' @return NULL (invisibly).
#' @export
print.s2index <- function(x, ...){
  cat("An s2index")
  return(invisible(NULL))
}

#' Look up IDs in a given s2index (beta version)
#'
#' Look up given s2cellids in a s2index and return the matching indices (without
#' replication). Expect changes to the API or even complete removal
#' of the function. Only for testing at the moment!
#'
#' @param id Cell ids in the s2cellid format.
#' @param index Index where the ids should be looked up.
#'
#' @return Vector of integers matching the s2cellids.
#' @export
s2lookup <- function(id, index){
  rslt <- NULL
  for(i in id){
    rslt <- c(rslt, index[[i]])
  }
  return(unique(rslt))
}
