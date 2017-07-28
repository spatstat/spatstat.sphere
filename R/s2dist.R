#' Distance on sphere
#'
#' @param x One or more points in a format acceptable to [globe::ensurelonlat()].
#' @param y One or more points in a format acceptable to [globe::ensurelonlat()].
#'
#' @details If both `x` and `y` contain more than one point they must contain
#' the same number of points and the distances between corresponding points are
#' returned.
#'
#' @return Numeric vector of distances.
#' @export
s2dist <- function(x, y){
  x <- globe::ensurelonlat(x)
  x$lat <- x$lat/180*pi
  x$lon <- x$lon/180*pi
  nx <- length(x$lon)
  y <- globe::ensurelonlat(y)
  y$lat <- y$lat/180*pi
  y$lon <- y$lon/180*pi
  ny <- length(y$lon)
  if(nx == 1L){
    x$lat <- rep(x$lat, ny)
    x$lon <- rep(x$lon, ny)
    nx <- ny
  }
  if(ny == 1L){
    y$lat <- rep(y$lat, nx)
    y$lon <- rep(y$lon, nx)
    ny <- nx
  }
  if(nx != ny){
    stop("When more than 1 location is supplied for both x and y they should have equal length.")
  }
  s2distHaversine(x$lat, x$lon, y$lat, y$lon)
}

s2distHaversine <- function(lat1, long1, lat2, long2){
  2 * asin(sqrt(sin((lat1 - lat2)/2)^2 + cos(lat1)*cos(lat2)*sin((long1 - long2)/2)^2))
}

pairdist.s2pp <- function(X, ...){
  x <- as.data.frame(X, format = "lon,lat")
  n <- nrow(x)
  lat <- x$lat/180*pi
  lon <- x$lon/180*pi
  distfun <- function(i,j){
    s2distHaversine(lat[i], lon[i], lat[j], lon[j])
  }
  rslt <- outer(1:n, 1:n, distfun)
  rad <- s2radius(X)
  if(rad != 1)
    rslt <- rslt * rad
  return(rslt)
}

closepairs.s2pp <- function(X, rmax, ..., chordal = FALSE){
  XX <- s2coords(X)
  rad <- s2radius(X)
  rmax <- rmax/rad
  if(!chordal)
    rmax <- 2*sin(rmax/2)
  closelist <- closepairs.pp3(pp3(XX[,1], XX[,2], XX[,3], box3(c(-1,1))), rmax = rmax, ...)
  if(!chordal)
    closelist$d <- asin(closelist$d/2)
  closelist$d <- 2*rad*closelist$d
  return(closelist)
}

nndist.s2pp <- function(X, ..., k = 1, chordal = FALSE){
  XX <- s2coords(X)
  rad <- s2radius(X)
  d <- nndist.pp3(pp3(XX[,1], XX[,2], XX[,3], box3(c(-1,1))), ..., k = k)
  if(!chordal)
    d <- asin(d/2)
  d <- 2*rad*d
  return(d)
}
