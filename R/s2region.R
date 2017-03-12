#' Generate a sphere
#'
#' @param radius Numeric. Defaults to 1. For convenience this can also be an
#' object of class `"s2region"` from where the radius and unitname then will be
#' extracted.
#' @param unitname Optional. Name of unit of length. Either a single character
#'   string, or a vector of two character strings giving the singular and plural
#'   forms, respectively.
#' @param ... Ignored.
#'
#' @return Object of class \code{"s2"} and super class \code{"s2region"}.
#' @export
s2 <- function(radius = 1, unitname = NULL, ...){
  if(inherits(radius, "s2region")){
    if(is.null(unitname))
      unitname <- radius$units
    radius <- s2radius(radius)
  }
  unitname <- as.units(unitname)
  s <- list(radius = radius, units = unitname)
  class(s) <- c("s2", "s2region")
  return(s)
}

#' Sphere approximating earth
#'
#' @return sphere of class \code{"s2"}.
#' @export
#'
#' @examples
#' earth <- s2earth()
#' area(earth)
s2earth <- function(){
  s2(radius = 6378.1, unitname = "km")
}

#' Radius of sphere
#'
#' @param x Object of class \code{"s2region"}.
#' @param format String to return the result as a numeric (`"numeric"`), a string
#'   with unitname (`"string"`) or a list with entries `radius` and `units`
#'   (`"list"`)
#'
#' @return numeric (default), string or list with the radius of the
#'   underlying sphere of the `"s2region"` and for string/list also the unitname.
#' @export
s2radius <- function(x, format = c("numeric", "string", "list")){
  format <- match.arg(format)
  rslt <- 1
  if(inherits(x, "s2pp"))
    x <- s2region(x)
  if(is.list(x) && !is.null(x$radius))
    rslt <- x$radius
  if(format == "string")
    rslt <- paste(rslt, ifelse(rslt==1, x$units$singular, x$units$plural))
  if(format == "list")
    rslt <- list(radius = rslt, units = x$units)
  return(rslt)
}

#' Extract spherical region
#'
#' @param x Object of class \code{"s2pp"}.
#'
#' @return object of class `"s2region"`.
#' @export
s2region <- function(x){
  verifyclass(x, "s2pp")
  return(x$domain)
}

#' Print object of class s2region
#'
#' @param x Object of class \code{"s2region"}.
#' @param ... Ignored.
#'
#' @return NULL (invisibly)
#' @export
print.s2region <- function(x, ...){
  msg1 <- ""
  if(inherits(x, "s2"))
    msg2 <- "The full sphere"
  if(inherits(x, "s2cap"))
    msg2 <- "A cap on the sphere"
  if(inherits(x, "s2polygon"))
    msg2 <- "A polygon on the sphere"
  if(inherits(x, "s2lonlatbox"))
    msg2 <- "A lon,lat box on the sphere"
  msg3 <- paste("of radius", s2radius(x, format = "string"))
  splat(msg1, msg2, msg3)
  return(invisible(NULL))
}

#' Test whether points are inside s2region
#'
#' @inheritParams s2pp
#' @param approx Logical to use approximate containment, which allows points
#' very slightly outside a polygonal region, but also ensures that points on the
#' border of a polygon are always accepted. Does not affect spherical caps at
#' the moment.
#'
#' @return Logical vector of same length as the number of rows in \code{coords}
#' @export
s2contains <- function(coords, region, approx = TRUE){
  if(!inherits(region, "s2region"))
    stop("Argument ", sQuote("region"), " must be a ", sQuote("s2region"), " object.")
  coords <- make_s2coords(coords)
  if(inherits(region, "s2polygon"))
    return(s2::S2Polygon_Contains(region, coords, approx = approx))
  if(inherits(region, "s2cap"))
    return(s2::S2Cap_Contains(region, coords))
  if(inherits(region, "s2"))
    return(rep(TRUE, nrow(coords)))
  stop("Can't check point containment for this type of s2region.")
}

#' Define spherical cap
#'
#' @param axis Axis (or centre) of the cap in a format acceptable to [`s2point`].
#' @param height Numeric of length one defining the height of the cap. Should be
#'   positive and at most twice the radius of the sphere. Not compatible with
#'   \code{dist}.
#' @param dist Numeric of length one defining the extend of the cap. All points
#'   on the sphere within distance \code{dist} of the point defined by
#'   \code{axis} are contained in the cap. Not compatible with \code{height}.
#' @inheritDotParams s2
#' @param simplify Logical. Whether to simplify to full sphere if the cap
#'   extends over the full sphere.
#'
#' @return Object of class `"s2cap"` and \code{"s2region"}.
#' @export
s2cap <- function(axis, height, dist, ..., simplify = TRUE){
  if(!missing(height) & !missing(dist))
    stop("Incompatible arguments ", sQuote("height"), " and ", sQuote("dist"),
         " both supplied.")
  # Resolve the containing sphere
  sphere <- do.call(s2, list(...))
  rad <- s2radius(sphere)
  # Resolve height missing or NULL
  if(missing(height) || is.null(height)){
    # Get info from dist
    if(missing(dist) || is.null(dist)){
      warning("No arguments specifying extend of cap. Returning an empty cap.")
      height <- -1
    } else{
      stopifnot(is.numeric(dist) && length(dist) == 1 && dist >= 0)
      dist <- dist/rad
      # Check dist is at most circumfrence/2
      if(dist >= pi)
        warning("dist is bigger than the maximal distance on the given sphere. Returning a full sphere.")
      # Convert to height
      height <- 1-cos(dist)
    }
  } else{
    stopifnot(is.numeric(height) && length(height) == 1 && height >= 0)
    height <- height/rad
    if(height >= 2 && simplify)
      warning("height is bigger than diameter of the given sphere. Returning a full sphere.")
  }
  # Return sphere if simplify is TRUE
  if(height >= 2 && simplify)
    return(sphere)
  axis <- s2point(axis)
  # Wrap up
  cap <- list(height = height, axis = axis, radius = rad, units = sphere$units)
  class(cap) <- c("s2cap", "s2region")
  return(cap)
}

#' Area of s2region
#'
#' @param w Object of class \code{"s2region"}.
#'
#' @examples
#' region <- s2(radius = 2)
#' area(region)
#' region <- s2cap(axis = c(lon=0,lat=0), height = 1, radius = 2)
#' area(region)
#' @export
area.s2region <- function(w){
  R <- w$radius
  if(inherits(w, "s2"))
    return(4*pi*R^2)
  if(inherits(w, "s2cap"))
    return(2*pi*R^2*max(0, w$height))
  if(inherits(w, "s2lonlatbox")){
    w <- list(lat = w$lat, lng = w$lon)
    return(R^2*s2::S2LatLngRect_area(w))
  }
  if(inherits(w, "s2polygon"))
    return(R^2*sum(w$areas * ifelse(w$holes, -1, 1)))
  stop("Unknown domain type.")
}

#' Create a s2polygon
#'
#' Create a s2polygon (object of class `"s2polygon"`).
#'
#' @param x List of loops (also called rings) to join as a polygon.
#' @param options List of options to control how the polygon is assembled. See
#' [`s2::S2Polygon`] for available options.
#' @inheritDotParams s2
#'
#' @export
#'
s2polygon <- function(x, ..., options = NULL){
  # Get options
  given <- names(options)
  valid <- names(formals(s2::S2Polygon))[-1]
  ok <- is.element(given, valid)
  if(any(!ok))
    warning("Some given options are not valid arguments of s2::S2Polygon and will be ignored.")
  options <- options[ok]
  # Make sure x is a list of things
  if(is.data.frame(x) | is.matrix(x)){
    x <- list(x)
  }
  x <- lapply(x, make_s2coords)
  n <- sapply(x, nrow)
  bad <- (n < 3)
  if(any(bad)){
    if(all(bad))
      stop("Loops must have three or more points.")
    warning(paste(sum(bad), "of the provided loops discarded due to having less than three points."))
    x <- x[!bad]
  }
  ## Initiate result
  rslt <- s2(...)
  ## Construct polygon using s2 package
  poly <- do.call(s2::S2Polygon, c(list(x=x), options))
  ## Fix bound format
  names(poly$bound)[names(poly$bound) == "lng"] <- "lon"
  poly$bound <- poly$bound[c("lon", "lat")]
  rslt <- c(rslt, poly)
  class(rslt) <- c("s2polygon", "s2region")
  return(rslt)
}

#' Intersection or union of two spherical regions
#'
#' Intersection or union of two objects of class `"s2region"`.
#'
#' @param x,y Region on the sphere of class `"s2region"`.
#'
#' @return Region on the sphere of class `"s2region"`.
#' @aliases s2union
#' @export
#' @examples
#' loop1 <- cbind(lon = c(0,60,60,0), lat = c(-40,-40,40,40))
#' loop2 <- cbind(lon = c(30,90,90,30), lat = c(-40,-40,40,40))
#' poly1 <- s2polygon(loop1)
#' poly2 <- s2polygon(list(loop2))
#' poly3 <- s2intersect(poly1, poly2)
#' poly4 <- s2union(poly1, poly2)
#' plot(poly1, eye = c(lon=45, lat=0), col = "red", eps = pi/100)
#' plot(poly2, col = "green", add = TRUE, eps = pi/100)
#' plot(poly3, col = "blue", add = TRUE, eps = pi/100)
#' plot(poly4, col = "cyan", eps = pi/100)
s2intersect <- function(x, y){
  ## Check args
  if(!inherits(x, "s2region"))
    stop("Argument ", sQuote("x"), " must be a ", sQuote("s2region"), " object.")
  if(!inherits(y, "s2region"))
    stop("Argument ", sQuote("y"), " must be a ", sQuote("s2region"), " object.")
  if(!all.equal(s2radius(x), s2radius(y)) | !identical(x$units, y$units))
    stop("Regions are defined on spheres with different radii.")

  ## If one argument is a full sphere we can just return the other
  if(inherits(x, "s2") | (inherits(x, "s2cap") && x$height >= 2))
    return(y)
  if(inherits(y, "s2") | (inherits(y, "s2cap") && y$height >= 2))
    return(x)

  ## Only polygon-polygon intersection allowed at the moment.
  if(!inherits(x, "s2polygon") || !inherits(y, "s2polygon"))
    stop("Can only do intersections of two polygons at the moment.")

  rslt <- s2::S2Polygon_intersection(x$loops, y$loops)
  return(wraps2polygon(rslt, s2(x)))
}

#' @export
s2union <- function(x, y){
  ## Check args
  if(!inherits(x, "s2region"))
    stop("Argument ", sQuote("x"), " must be a ", sQuote("s2region"), " object.")
  if(!inherits(y, "s2region"))
    stop("Argument ", sQuote("y"), " must be a ", sQuote("s2region"), " object.")
  if(!isTRUE(all.equal(s2radius(x), s2radius(y))) | !identical(x$units, y$units))
    stop("Regions are defined on spheres with different radii or units.")

  ## If one argument is a full sphere we can just return it
  if(inherits(x, "s2") | (inherits(x, "s2cap") && x$height >= 2))
    return(x)
  if(inherits(y, "s2") | (inherits(y, "s2cap") && y$height >= 2))
    return(y)

  ## Identical regions can be returned
  if(identical(x, y))
    return(x)

  ## Only polygon-polygon union allowed at the moment.
  if(!inherits(x, "s2polygon") || !inherits(y, "s2polygon"))
    stop("Can only do unions of two polygons at the moment.")

  rslt <- s2::S2Polygon_union(x$loops, y$loops)
  return(wraps2polygon(rslt, s2(x)))
}

wraps2polygon <- function(poly, sphere){
  class(poly) <- c("s2polygon", "s2region")
  poly$radius <- s2radius(sphere)
  poly$units <- sphere$units
  return(poly)
}

#' Create a list of loops on the sphere
#'
#' Create a list of loops on the sphere.
#'
#' @param ... Input to create the loops from. Each input should be a matrix or
#' data.frame with the coordinates of each vertex of the loop in each row. A
#' `s2polygon` or `s2looplist`.
#'
#' @return A list of loops with class `s2looplist`.
#' @export
#'
#' @examples
#' loop1 <- cbind(lon = c(0,60,60,0), lat = c(-40,-40,40,40))
#' loop2 <- cbind(lon = c(30,90,90,30), lat = c(-40,-40,40,40))
#' x <- s2looplist(loop1, loop2)
#' y <- s2looplist(x, loop1)
#' z <- s2looplist(x, s2polygon(loop1))
s2looplist <- function(...){
  dots <- list(...)
  extractloops <- function(x){
    if(inherits(x, "s2cellid"))
      x <- s2cell(x)
    if(inherits(x, "s2cell"))
      return(x$vertices)
    if(inherits(x, "s2looplist"))
      return(x)
    if(inherits(x, "s2polygon"))
      return(x$loops)
    return(list(globe::ensure3d(x, single = FALSE)))
  }
  loops <- lapply(dots, extractloops)
  loops <- unlist(loops, recursive = FALSE)
  class(loops) <- "s2looplist"
  return(loops)
}

#' Generic function to extract or construct a box in longitude and latitude coordinates
#'
#' Generic function to extract or construct a box in longitude and latitude coordinates.
#'
#' @param ... Parameters passed to methods. For `s2lonlatbox.s2region` only
#'   first argument (named `x`) is used and the rest ignored. For
#'   `s2lonlatbox.default` the first two arguments (named `lon` and `lat`) are
#'   used to define the region and the rest are passed to [s2()] to define the
#'   radius and units of the underlying sphere.
#'
#' @return Box in longitude and latitude coordinates (of class `"s2lonlatbox"`).
#' @seealso s2lonlatbox.s2region s2lonlatbox.default
#' @export
#'
#' @examples
#' cap <- s2cap(c(1, 0, 0), h = 0.1)
#' s2lonlatbox(cap)
#' s2lonlatbox(lon = c(-30, 30), lat = c(45, 90))
#' ## Very large box crossing the 180th meridian
#' big <- s2lonlatbox(lon = c(1, -1), lat = c(-90, 90))
#' area(big)
#' ## Very small box **not** crossing the 180th meridian
#' small <- s2lonlatbox(lon = c(-1, 1), lat = c(-90, 90))
#' area(small)
#'
s2lonlatbox <- function(...){
  UseMethod("s2lonlatbox")
}

#' @describeIn s2lonlatbox Extraction method for `"s2region"`.
#' @param x Object of class `"s2region"`.
#' @export
s2lonlatbox.s2region <- function(x, ...){
  if(inherits(x, "s2lonlatbox"))
    return(x)
  rslt <- s2radius(x, format = "list")
  if(inherits(x, "s2")){
    rslt$lon <- c(-180, 180)
    rslt$lat <- c(-90, 90)
  } else if(inherits(x, "s2cap")){
    cap <- s2::S2Cap_GetRectBound(x)
    rslt$lon <- cap$lng
    rslt$lat <- cap$lat
  } else if(inherits(x, "s2polygon")){
    rslt$lon <- x$bound$lon
    rslt$lat <- x$bound$lat
  } else{
    stop("Can't extract s2lonlatbox for this s2region.")
  }
  class(rslt) <- c("s2lonlatbox", "s2region")
  return(rslt)
}

#' @describeIn s2lonlatbox Default method to generate `"s2lonlatbox"` from a
#'   longitude and latitude interval.
#' @param lon Numeric of length 2 specifying longitude interval in degrees. For
#'   boxes not crossing the 180 degree meridian the input should be such that
#'   `-180 <= lon[1] < lon[2] <= 180`. For boxes crossing the 180 degree meridian
#'   the input shoud be such that `-180 <= lon[2] < lon[1] <= 180`.
#' @param lat Numeric of length 2 with `-90 <= lat[1] < lat[2] <= 90` specifying
#'   latitude interval in degrees.
#' @export
s2lonlatbox.default <- function(lon = c(-180, 180), lat = c(-90, 90), ...){
  rslt <- s2(...)
  x <- s2::S2LatLngRect(c(lat[1], lon[1]), c(lat[2], lon[2]))
  rslt$lat <- x$lat
  rslt$lon <- x$lng
  class(rslt) <- c("s2lonlatbox", "s2region")
  return(rslt)
}
