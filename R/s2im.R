#' Make a grid over spherical polygonal region
#'
#' @param x Polygon on the sphere of class `s2polygon`.
#' @param level Integer level of the pixels between 0 (each pixel is huge --
#' millions of square kilometres on Earth) and 30 (each pixel is tiny --
#' approximately a square cm on Earth).
#' @param ... Ignored.
#'
#' @return Object of class "s2grid", which is just a vector of unique ids of
#' the cells making up the grid.
#' @export
#'
#' @examples
#' poly <- s2polygon(cbind(lon = c(-20,-10,10,20,10,-10), lat = c(0, -10, -10, 0, 10, 10)))
#' grid <- s2grid.s2polygon(poly, 5)
#' length(grid)
s2grid.s2polygon <- function(x, level, ...){
  grid <- s2::s2_covering_cell_ids(as.s2_polygon(x), min_level = level, max_level = level)
  grid <- unclass(grid)[[1]]
  class(grid) <- c("s2grid", class(grid))
  return(grid)
}

#' Make an Image/Raster Object on the Sphere
#'
#' @param grid An `s2grid` defining the cells/pixels of the image.
#' @param values A `numeric` of same length as `grid` containing the image
#' value at each cell.
#'
#' @return Object of class `s2im`.
#' @export
#'
#' @examples
#' poly <- s2polygon(cbind(lon = c(-20,-10,10,20,10,-10), lat = c(0, -10, -10, 0, 10, 10)))
#' grid <- s2grid.s2polygon(poly, 5)
#' im <- s2im(grid, 1:length(grid))
#' plot(im)
s2im <- function(grid, values){
  d <- data.frame(id = grid, values = values)
  class(d) <- "s2im"
  return(d)
}

#' Convert Function on the Sphere to an Image/Raster on the Sphere
#'
#' @param f function of degrees longitude and latitude in that order.
#' @param region Object of class `s2region` where the function should be evaluated.
#' @param ... Ignored.
#' @param level Integer level of the pixels between 0 (each pixel is huge --
#' millions of square kilometres on Earth) and 30 (each pixel is tiny --
#' approximately a square cm on Earth).
#' @return Object of class `s2im`.
#' @export
#'
#' @examples
#' poly <- s2polygon(cbind(lon = c(-20,-10,10,20,10,-10), lat = c(0, -10, -10, 0, 10, 10)))
#' f <- function(lon, lat){ sin(lon) }
#' im <- as.s2im.function(f, poly, level = 6)
#' plot(im)
as.s2im.function <- function(f, region = NULL, ..., level = NULL){
  grid <- s2grid.s2polygon(region, level)
  centers <- as.s2pp(grid)
  co <- coords(centers)
  co <- globe::ensurelonlat(co)
  val <- f(co$lon, co$lat)
  s2im(grid, val)
}

#' Convert a Grid on the Sphere to a Point Pattern on the Sphere
#'
#' @param X Object of class `s2grid`.
#' @param ... Ignored.
#'
#' @return Point pattern of class `s2pp`. The points are the grid cell
#' centers.
#' @export
#'
#' @examples
#' poly <- s2polygon(cbind(lon = c(-20,-10,10,20,10,-10), lat = c(0, -10, -10, 0, 10, 10)))
#' grid <- s2grid.s2polygon(poly, 5)
#' pp <- as.s2pp(grid)
as.s2pp.s2grid <- function(X, ...){
  centers <- s2::s2_cell_center(X)
  co <- data.frame(lon = s2::s2_x(centers), lat = s2::s2_y(centers))
  s2pp(co)
}

#' Convert an Image/Raster on the Sphere to a Marked Point Pattern on the Sphere
#'
#' @param X Object of class `s2im`.
#' @param ... Ignored.
#'
#' @return Marked point pattern of class `s2pp`. The points are the image cell
#' centers and the marks are the corresponding values.
#' @export
#'
#' @examples
#' poly <- s2polygon(cbind(lon = c(-20,-10,10,20,10,-10), lat = c(0, -10, -10, 0, 10, 10)))
#' f <- function(lon, lat){ sin(lon) }
#' im <- as.s2im.function(f, poly, level = 6)
#' pp <- as.s2pp(im)
#' plot(pp)
as.s2pp.s2im <- function(X, ...){
  centers <- s2::s2_cell_center(X$id)
  co <- data.frame(lon = s2::s2_x(centers), lat = s2::s2_y(centers))
  s2pp(co, marks = X$values)
}

#' Plot Image/Raster on the Sphere
#'
#' @param x Object of class `s2im`.
#' @param ... passed to <plot.s2pp>.
#' @export
#'
#' @examples
#' poly <- s2polygon(cbind(lon = c(-20,-10,10,20,10,-10), lat = c(0, -10, -10, 0, 10, 10)))
#' f <- function(lon, lat){ sin(lon) }
#' im <- as.s2im.function(f, poly, level = 6)
#' plot(im)
plot.s2im <- function(x, ...){
  X <- as.s2pp.s2im(x)
  plot(X, ...)
}

#' Print Method for Image/Raster on the Sphere
#'
#' @param x Object of class `s2pp`
#' @param ... Ignored.
#'
#' @return NULL
#' @export
print.s2im <- function(x, ...){
  splat("An s2im.")
}
