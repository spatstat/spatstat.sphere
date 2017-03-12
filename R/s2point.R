#' @title Single point on the (unit) sphere
#'
#' @description Make an object of class `"s2point"` representing a single point
#' on the unit sphere.
#'
#' @param ... Input to interpret as a point on the sphere. Either in terms of
#' three Eucledian coordinates (x,y,z) or two spherical coordinates
#' (latitude,longitude). See Details.
#'
#' @details For three Eucledian coordinates the following type of inputs are
#' recognised:
#' - Three inputs with a single real number (optionally named `x`, `y` and `z`).
#' - A list or data.frame with three real numbers as elements (optionally named `x`, `y` and `z`).
#' - A single vector or matrix with three real numbers.
#'
#' For two spherical coordinates a named input is required. The names can either
#' be given as argument names of the function or as the names of the provided
#' input. Input latitude is partially matched to `"latitude"` and input longitude
#' is either partially matched to `longitude` or exactly matched to `lng`. In all
#' cases measurements are in degrees and should be in \eqn{[-90,90]} (latitude)
#' and \eqn{[-180,180]} (longitude).
#' The following type of inputs are recognised:
#' - Two named inputs with a single real number.
#' - A list or data.frame with two named elements containing a real number each.
#' - A single vector with two named elements containing a real number each.
#' - A two column matrix with a single row of real numbers and a `colnames`
#' attribute.
#'
#' @return Object of class \code{"s2pp"}.
#' @export
#'
#' @examples
#' longitude <- c(90)
#' latitude <- c(45)
#' xyz <- s2point(lon = longitude, lat = latitude)
#' x <- xyz[1]
#' y <- xyz[2]
#' z <- xyz[3]
#' xyz2 <- s2point(x=x, y=y, z=z)
#'
s2point <- function(...){
  input <- list(...)
  n <- length(input)
  if(n<1 | n>3)
    stop("Please provide one, two or three inputs.")
  if(length(input) == 1)
    input <- input[[1]]
  # Now input should be length 2 or 3
  n <- length(input)
  if(n < 2 | n > 3)
    stop("Could not interpret output as containing two or three coordinates.")
  # Initialize output object
  output <- NULL
  if(is.matrix(input)){
    if(nrow(input) > 1){
      input <- t(input)
    }
    input <- as.data.frame(input)
  }
  ## Eucledian input
  if(n == 3){
    input <- unlist(input)
    nam <- tolower(names(input))
    if(setequal(nam, c("x", "y", "z"))){
      names(input) <- nam
      input <- input[c("x", "y", "z")]
    }
    norm <- sqrt(sum(input^2))
    if(norm<.Machine$double.eps)
      stop("Provided Eucledian coordinates almost all zeros. Can't normalize unit vector.")
    output <- input/norm
  }
  ## Spherical input
  if(n == 2){
    input <- unlist(input)
    if(is.null(names(input)))
      stop("Named arguments must be provided for latitude and longitude coordinates.")
    names(input) <- tolower(names(input))
    nam <- names(input)
    # Recover lng if given
    nam[nam=="lng"] <- "lon"
    lonlatorder <- pmatch(nam, c("longitude", "latitude"))
    if(!setequal(lonlatorder, 1:2)){
      stop("Couldn't resolve latitude and longitude from names. Check spelling?")
    }
    input <- input[lonlatorder]
    if(input[1] < -180 | input[1] > 180)
      stop("Longitude outside [-180,180].")
    if(input[2] < -90 | input[2] > 90)
      stop("Latitude outside [-90,90].")
    output <- globe::ensure3d(input)
  }
  if(!is.null(output)){
    names(output) <- c("x", "y", "z")
    class(output) <- "s2point"
    return(output)
  }
  stop("Couldn't interpret input.")
}

#' Print object of class s2point
#'
#' @param x Object of class `"s2point"`.
#' @param ... Ignored.
#' @return NULL (invisibly)
#' @export
print.s2point <- function(x, ...){
  splat("Point on the sphere represented by the unit direction vector:")
  print(unclass(x))
  splat(paste("Corresponding to spherical/geographical coordinates:"))
  print(unlist(globe::ensurelonlat(x)))
  return(invisible(NULL))
}
