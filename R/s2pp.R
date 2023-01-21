make_s2coords <- function(x){
  if(inherits(x, "s2pp")){
    return(coords(x))
  }
  if(inherits(x, "s2_geography")){
    x <- cbind(lon = s2::s2_x(x), lat = s2::s2_y(x))
  }
  lat_names <- c("lat", "latitude")
  lon_names <- c("lon", "long", "longitude", "lng")
  # Handle matrix
  if(is.matrix(x)){
    is_latlng <- any(colnames(x) %in% c(lat_names, lon_names))
    # Detect correct format and return directly
    if( ncol(x) == 3 && !is_latlng){
      colnames(x) <- c("x", "y", "z")
      return(x)
    }
    # Otherwise convert to data.frame to handle all in one place
    x <- as.data.frame(x)
  }
  # Handle list/data.frame:
  if(is.list(x)){
    nam <- names(x)
    is_latlng <- any(nam %in% c(lat_names, lon_names))
    # Assume x,y,z coordinates:
    if(length(x) == 3 && !is_latlng){
      if( length(x[[1]]) != length(x[[2]]) | length(x[[1]]) != length(x[[3]]) )
        stop("Unequal length of x,y,z coordinate entries.")
      return(cbind(x = x[[1]], y = x[[2]], z = x[[3]]))
    }
    # Assume lat,long coordinates
    if(length(x) == 2){
      lat <- x[[match.arg(nam, lat_names, several.ok = TRUE)[1]]]
      lon <- x[[match.arg(nam, lon_names, several.ok = TRUE)[1]]]
      if(is.null(lat) | is.null(lon)){
        stop("Names of input couldn't be matched to latitude and longitude. Check spelling?")
      }
      x <- globe::ensure3d(data.frame(lon, lat), single = FALSE)
      colnames(x) <- c("x", "y", "z")
      return(x)
    }
  }
  stop("Input must be a matrix or data.frame with either 2 or 3 columns.")
}

#' @title Spherical point pattern
#'
#' @description  Make an object of class `"s2pp"` representing a point pattern
#' over a region on the sphere (of class `"s2region"`).
#'
#' @param coords `data.frame`, `matrix` or `list`. Object containing spherical
#'   coordinates of the points. This must be either longitude and latitude in
#'   degrees or 3D Cartesian coordinates. In the former case two named columns
#'   or list entries are required. In the latter three colums or list entries
#'   are required and assumed to be in x, y, z order (disregarding any names).
#' @param region Region on the sphere where the points occur. Object of class
#'   `"s2region"`.
#' @param marks Marks attached to the points.
### #' @param  ... Arguments passed to [`s2`] if the region is not otherwise specified.
#' @inheritDotParams s2
#' @param check Logical. Check that points actually are inside the provided
#'   `region`.
#'
#' @return Object of class \code{"s2pp"}.
#' @export
#'
#' @examples
#' long <- c(90, 180, -90)
#' lat <- c(45, 0, -45)
#' coords <- data.frame(long, lat)
#' X <- s2pp(coords)
#'
s2pp <- function(coords, region = NULL, marks = NULL, ..., check = TRUE){
  if(missing(coords)){
    coords <- matrix(numeric(0), ncol=3)
  }
  if(is.null(region)){
    region <- do.call(s2, list(...))
  }
  coords <- make_s2coords(coords)
  ## Coordinate types:
  ctype <- rep("s", 3)
  n <- nrow(coords)
  if (check && n > 0) {
    ok <- s2contains(coords, region)
    nout <- sum(!ok)
    if (nout > 0) {
      warning(paste(nout, ngettext(nout, "point was", "points were"),
                    "rejected as lying outside the specified window"))
      rejects <- coords[!ok,]
      coords <- coords[ok,]
      n <- nrow(coords)
    }
  }
  else{
    nout <- 0
  }
  if(!is.null(marks)){
    if(!is.data.frame(marks)){
      marks <- data.frame(marks = marks)
    }
    if(nout>0){
      marks <- marks[ok,]
    }
    ctype <- c(ctype, rep("mark", ncol(marks)))
    coords <- cbind(coords, marks)
  }
  rslt <- ppx(coords, domain = region, coord.type = ctype)
  class(rslt) <- c("s2pp", class(rslt))

  if (check && anyDuplicated(rslt))
    warning("data contain duplicated points")
  if (nout > 0)
    attr(rslt, "rejects") <- rejects
  return(rslt)
}

#' Print object of class s2pp
#'
#' @param x Object of class `"s2pp"`.
#' @param ... Ignored.
#' @return NULL (invisibly)
#' @export
print.s2pp <- function(x, ...){
  splat(paste("Spherical point pattern with", npoints(x), "points."))
  splat("Region (of class 's2region'):")
  print(s2region(x))
  # Temporary hack: Detect attribute set by simulation algorithm for DPPs
  nmean <- attr(x, "nmean")
  if(!is.null(nmean)){
    cat(paste("It has been simulated from a model where\n",
              "the expected number of points is: ",
              signif(nmean,4), "\n", sep=""))
  }
  return(invisible(NULL))
}

#' Extract coordinates of a point pattern on a sphere
#'
#' @param x Object of class \code{"s2pp"}
#' @param ... Ignored
#'
#' @return a \code{data.frame} with one row for each point, containing the
#'   coordinates.
#' @export
#'
coords.s2pp <- function(x, ...){
  return(s2coords(x))
}

#' Extract Coordinates of a Point Pattern on a Sphere
#'
#' @param x Object of class \code{"s2pp"}
#'
#' @return a \code{data.frame} with one row for each point, containing the
#'   coordinates.
#' @export
#'
s2coords <- function(x){
  verifyclass(x, "s2pp")
  return(coords.ppx(x, temporal = FALSE, local = FALSE))
}

#' Export coordinates of a spherical point pattern into a data.frame
#'
#' @param x Spherical point patten of class \code{"s2pp"}.
#' @param ... Ignored.
#' @param format String. Type of coordinates used for return coordinates.
#'   Defaults to `"long,lat"`.
#'
#' @details Return format would usually be a latitude,longitude format (in any
#' order) by specifying the desired names of the two output columns. Valid
#' choices for latitude are `lat` and `latitude`. Valid choices for longitude
#' are `lon`, `long`, `lng` and `longitude`. Alternatively, any order of x,y,z
#' coordinates are obtained by supplying a commaseparated format like `"y,z,x"`.
#'
#' @return \code{data.frame} with coordinates.
#' @export
#'
as.data.frame.s2pp <- function(x, ..., format = "long,lat"){
  x <- s2coords(x)
  format <- strsplit(format, ",", fixed = TRUE)[[1]]
  if(setequal(format, c("x", "y", "z"))){
    return(x[format])
  }
  lat_names <- c("lat", "latitude")
  lon_names <- c("lon", "long", "longitude", "lng")
  ## Check length 2 format with one lat and one lon
  if(length(format)!=2 || !any(format %in% lat_names) || !any(format %in% lon_names))
    stop("Wrong format string.")
  x <- globe::ensurelonlat(x)
  ## Switch order if lat is first
  if(format[1] %in% lat_names)
    x <- x[2:1]
  x <- as.data.frame(x)
  names(x) <- format
  return(x)
}

#' Convert data to class s2pp
#'
#' Tries to coerce any reasonable kind of data to a spherical point pattern (an
#' object of class \code{"s2pp"}) for use by the \pkg{spatstat.sphere} package).
#'
#' Converts the dataset \code{X} to a point pattern (an object of class
#' \code{"s2pp"}).
#'
#'
#### Unused old documention bit #####
# This function is normally used to convert an existing point pattern dataset,
# stored in another format, to the \code{"s2pp"} format.  To create a new point
# pattern from raw data such as lat,long coordinates, it is normally easier to
# use the creator function \code{\link{s2pp}}.
#
# The dataset \code{X} may be:
#
# \itemize{
#   \item an object of class \code{"s2pp"}
#   \item a matrix or data frame with at least two columns
#   \item a structure with entries \code{lat}, \code{long} which are numeric
#     vectors of equal length
#   \item a numeric vector of length 2, interpreted as the coordinates of a
#     single point.
# }
#
# The first case is typically used to change (or ensure) a specific coordinate
# format by specifying either the argument \code{region} directly or indirectly
# using the argument \code{coordsys} which is passed to \code{\link{sphere}}
# through the additional arguments \dots. In the last three cases the default
# behaviour is to assume the region is the entire sphere with coordinates of
# type \code{"geo_deg"} (see \code{\link{s2pp}}). Alternatively the region can
# be specified by the argument \code{region} or through the additional
# arguments \dots.
#
# If \code{X} is a matrix or data frame, the first and second columns will be
# interpreted as the \eqn{lat} and \eqn{long} coordinates respectively. Any
# additional columns will be ignored.
####################################
#' The function \code{as.s2pp} is generic, with methods for the classes
#' \code{"s2pp"}, \code{"matrix"}, \code{"data.frame"} and a default method.
#'
#' Point pattern datasets can also be created by the function \code{\link{s2pp}}.
#'
#' @param X Data which will be converted into a point pattern
#' @inheritParams s2pp
#' @inheritDotParams s2
#' @return An object of class \code{"s2pp"} (see \code{\link{s2pp}}) describing
#'   the spherical point pattern and its region.
#' @keywords spatial manip
#' @export
as.s2pp <- function(X, ...){
  UseMethod("as.s2pp")
}

#' @describeIn as.s2pp Method for s2pp objects.
#' @export
as.s2pp.s2pp <- function(X, ...){
  verifyclass(X, "s2pp")
  return(X)
}

#' @describeIn as.s2pp Method for data frames.
#' @export
as.s2pp.data.frame <- function(X, ..., region = NULL){
  if(ncol(X)<2) stop("Need at least 2 columns for spherical coordinates.")
  if(is.null(region))
    region <- do.call(s2, list(...))
  co <- make_s2coords(X)
  s2pp(X, region = region, check = TRUE)
}

#' @describeIn as.s2pp Method for matrices.
#' @export
as.s2pp.matrix <- as.s2pp.data.frame

#' @describeIn as.s2pp Default method.
#' @export
as.s2pp.default <- function(X, ..., region = NULL){
  if(is.null(region))
    region <- do.call(s2, list(...))
  s2pp(s2coords(X), region = region, check = TRUE)
}

#' Extract a subset of a point pattern on a sphere.
#'
#' Extract a subset of a point pattern on a sphere. Extraction of a subset has
#' the effect of thinning the points.
#'
#' @param x object of class \code{"s2pp"}.
#' @param i Subset index. A valid subset index in the usual R sense, indicating
#'   which points should be retained.
#' @param j Ignored. (Required for compatibility with the generic function.)
#' @param drop Ignored. (Required for compatibility with the generic function.)
#' @param \dots Ignored. (Required for compatibility with the generic function.)
#' @param clip Logical value indicating how to form the s2region of the
#' resulting point pattern, when i is a window. If clip=FALSE (the default), the
#' result has s2region equal to i. If clip=TRUE, the resulting s2region is the
#' intersection between the s2region of x and the s2region i.
#' @return object of class \code{"s2pp"}.
#' @export
#' @examples
#'
#' X <- s2pp(data.frame(lat = c(-45,0,45), long = c(-10,0,160)))
#' X[1:2]
#'
"[.s2pp" <- function(x, i, j, drop, ..., clip = FALSE) {
  verifyclass(x, "s2pp")
  if(missing(i) || npoints(x) == 0)
    return(x)
  if(inherits(i, "s2region")){
    region <- i
    if(clip)
      region <- s2intersect(region, s2region(x))
    ok <- s2contains(s2coords(x), region)
    x <- s2pp(s2coords(x)[ok,], region = region, check = FALSE)
  } else{
    x$data <- x$data[i,]
  }
  return(x)
}

s2borderdist <- function(x, region = NULL){
  if(inherits(x, "s2pp")){
    region <- s2region(x)
    x <- s2coords(x)
  }
  if(is.null(region))
    stop("An s2region must be supplied.")
  x <- globe::ensure3d(x, single = FALSE)
  if(inherits(region, "s2polygon")){
    x <- s2::as_s2_point(x)
    ## Loops as lines
    loop_lines <- lapply(region$loops, function(x) s2::s2_make_line(x[,1], x[,2]))
    ## List giving the distance from each point to each loop
    loopdist <- lapply(loop_lines, s2::s2_distance, x = x, radius = s2radius(region))
    ## Shortest border distance for each point
    loopdist <- Reduce(pmin, loopdist)
    return(loopdist)
    # return(s2::S2Polygon_border_dist(x, region$loops))
  }
  if(inherits(region, "s2") | inherits(region, "s2cap") && region$height >= 2*s2radius(region))
    return(rep(Inf, nrow(x)))
  stop("Can only handle polygons and full spheres at the moment.")
}

#' Superimpose several point patterns on the sphere
#'
#' @param ... Any number of arguments, each of which represents a point pattern
#'   on the sphere.
#' @param region Optional. Object of class `"s2region"` determining the region
#'   for the resulting pattern.
#' @param check Logical value passed to [s2pp()] determining whether to check
#'   the geometrical validity of the resulting pattern.
#'
#' @return Point pattern on the sphere of class `"s2pp"`.
#'
#' @importFrom spatstat.geom superimpose
#' @exportS3Method superimpose s2pp
#' @export superimpose.s2pp
#'
#' @examples
#' p1 <- s2runif(10)
#' p2 <- s2runif(20)
#' s1 <- superimpose(p1, p2)
superimpose.s2pp <- function(..., region = NULL, check = TRUE){
  Xlist <- lapply(list(...), as.s2pp)
  coords <- lapply(Xlist, s2coords)
  coords <- Reduce(rbind, coords)
  if(is.null(region)){
    regions <- lapply(Xlist, s2region)
    region <- Reduce(s2union, regions)
    check <- FALSE
  }
  return(s2pp(coords, region, check = check))
}
