#' Make s2cellids from points on the sphere
#'
#' Create a vector of s2cellids corresponding to the cells at the given level
#' containing the given points. The default level (30) corresponds to leaf
#' cells (finest/maximal level).
#'
#' @param x Points on the sphere in a standard format.
#' @param level Integer between 0 and 30 (incl).
#'
#' @return A character vector with s2cellids with the additional class `s2cellid`.
#' @export
s2cellid <- function(x, level = 30){
  x <- globe::ensure3d(x, single = FALSE)
  x <- matrix(x, ncol = 3)
  id <- s2::S2CellIdFromPoint(x, level)
  class(id) <- "s2cellid"
  return(id)
}

#' Extract s2 cell level from data
#'
#' @param x object inheriting from `s2cellid`
#' @param simplify logical to return a single value if all levels are equal.
#'
#' @return Integer vector.
#' @export
#'
#' @examples
#' x <- rbind(c(1,0,0), c(0,1,0))
#' id <- s2cellid(x)
#' s2level(id)
#' s2level(id, simplify = TRUE)
#'
#' id <- s2cellid(x, c(1,30))
#' s2level(id)
s2level <- function(x, simplify = FALSE){
  id <- x$id
  n <- length(id)
  lev <- x$lev
  if(!simplify && n > 1L){
    if(length(lev) != n){
      lev <- rep(lev, n)
    }
  }
  if(simplify && n > 1L){
    lev <- unique(lev)
    if(length(lev) > 1L){
      stop("Different levels present -- can't simplify to one value!")
    }
  }
  return(lev)
}

#' Coerce multiple inputs to a s2cellid vector
#'
#' @param ... Inputs. Currently only objects of class `"s2cellid"` are understood.
#'
#' @return Object of class `"s2cellid"`.
#' @export
as.s2cellid <- function(...){
  dots <- list(...)
  ok <- sapply(dots, inherits, what = "s2cellid")
  if(!any(ok))
    stop("No valid input provided.")
  ids <- dots[ok]
  lev <- lapply(ids, s2level)
  ids <- lapply(ids, function(x) x$id)
  rslt <- list(id = unlist(ids), level = unlist(lev))
  class(rslt) <- "s2cellid"
  return(rslt)
}

#' Print object of class s2cellid
#'
#' @param x Object of class `"s2cellid"`.
#' @param ... Ignored.
#' @return NULL (invisibly)
#' @export
print.s2cellid <- function(x, ...){
  splat("Object of class s2cellid containing", length(x$id), "ids.")
}

#' Make a list of s2cells
#'
#' Make a list of s2cells
#'
#' @param x Input to create cells from. Currently only a vector of s2cellid
#' tokens are supported.
#'
#' @export
s2cell <- function(x){
  rslt <- s2::S2Cell(x)
  class(rslt$vertices) <- "s2looplist"
  rslt$centers <- s2::S2CellId_ToPoint(x)
  class(rslt) <- "s2cell"
  return(rslt)
}

#' Retreive string representation of a s2cell
#'
#' Retreive string representation of a s2cell
#'
#' @param x Input to retrieve string representations from.
#' Current options are a vector of s2cellid tokens, a list of s2cells, or a
#' three-column matrix representing points on the sphere.
#' @param level Integer between 0 and 30 (incl.) to specify the desired s2cell
#' level. Defaults to 30 for point data and is ignored for s2cell(id) data.
#' @param binary Logical to request binary representation where four subcells
#' are represented as 00, 01, 10 and 11 rather than by 0, 1, 2, 3 which is the
#' default.
#' @param zerotail Logical to turn on/off zero-padding at the end when `binary = TRUE`.
#'
#' @export
s2cellstring <- function(x, level = 30L, binary = FALSE, zerotail = binary){
  if(inherits(x, "s2cell"))
    x <- x$ids
  if(!inherits(x, "s2cellid")){
    ## Now we assume x is in point form
    x <- s2cellid(x, level = level)
  }
  x <- s2::S2CellId_ToString(x$id)
  if(!binary)
    return(x)
  str2bin <- function(y, width = 2){
    y <- sub("0", ifelse(width==2, "00", "000"), y)
    y <- sub("1", ifelse(width==2, "01", "001"), y)
    y <- sub("2", ifelse(width==2, "10", "010"), y)
    y <- sub("3", ifelse(width==2, "11", "011"), y)
    y <- sub("4", "100", y)
    y <- sub("5", "101", y)
    y
  }
  y <- strsplit(x, "")[[1]]
  face <- str2bin(y[1], 3)
  y <- y[-(1:2)]
  y <- paste(c(face,sapply(y, str2bin)), collapse = ",")
  y <- paste0(y, ",1")
  if(zerotail && level < 30L){
    y <- paste(y, ifelse(level == 29, "00", "0..0"), sep = ",")
  }
  y
}

## #' Plot the Outline of s2cells on the sphere
## #'
## #' Plot the outline of s2cells of class `"s2cell"` in base
## #' graphics via the [globe package][globe::globe-package].
## #'
## #' @param x Object of class `"s2cell"` to plot.
## #' @param ... parameters passed to [plot.s2polygon()], which does the actual
## #' plotting.
## #' @param radius Positive real number. Optional radius to use for the sphere
## #' (only effects how to interpret the argument `eps` in [plot.s2polygon()]).
## #'
## #' @return NULL (invisibly)
## #' @export
## plot.s2cell <- function(x, ...){
##   nam <- names(x)
##   nam[nam=="vertices"] <- "loops"
##   names(x) <- nam
##   x$radius <- 1
##   plot.s2polygon(x, ...)
## }

#' Approximate a region on the sphere by a covering of s2cells
#'
#' Approximate a region on the sphere by a (possibly interior) covering of
#' s2cells.
#'
#' @param x Region to cover. Currently it must be a polygon, cap or full sphere.
#' @param max_cells Positive integer. Maximal number of cells to use in the
#' covering.
#' @param min_level Integer between 0 and 30 specifying the lowest cell level to
#' use. Must be less than or equal to `max_level`.
#' @param max_level Integer between 0 and 30 specifying the highest cell level to
#' use. Must be greater than or equal to `min_level`.
#' @param interior Logical to get an interior covering.
#' @return A vector of s2cellids (of class `s2cellid`).
#' @export
#' @examples
#' # Covering of entire sphere at level 1
#' s2covering(s2earth(), min_level = 1, max_level = 1)
#'
s2covering <- function(x, max_cells = 8, min_level = 0, max_level = 30, interior = FALSE){
  if(inherits(x, "s2"))
    x <- s2cap(c(1,0,0), height = 2, simplify = FALSE)
  class(x) <- class(x)[1]
  rslt <- s2::S2Covering(x = x, max_cells = max_cells, min_level = min_level,
                         max_level = max_level, interior = interior)
  class(rslt) <- "s2cellid"
  return(rslt)
}
