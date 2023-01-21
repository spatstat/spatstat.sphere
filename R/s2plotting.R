#' Plot a spherical point pattern
#'
#' Plot a spherical point pattern of class `"s2pp"` in base graphics via the
#' [globe package][globe::globe-package].
#'
#' @param x Point pattern of class `"s2pp"` to plot.
#' @param ... parameters passed to [globe::globepoints()] to control view point
#' (`eye`), orientation (`top`), and grapical parameters such as point color
#' (`col`), point size (`cex`), point type (`pch`) etc.
#' @param add Logical to add the points to existing plot.
#' @param region Either logical to add the boundary of the `"s2region"` of the
#' point pattern or a user supplied `"s2region"` to which the point pattern will
#' be restricted.
#' @param longrid Numeric defining a grid of longitude lines to be plotted if
#' `add = FALSE`. Either a vector of numbers between -180 and 180 or a single
#' numeric used as a stepsize between consecutive longitude lines. Value of zero
#' or `NULL` disables longitude lines.
#' @param latgrid Numeric defining a grid of latitude lines to be plotted if
#' `add = FALSE`. Either a vector of numbers between -90 and 90 or a single
#' numeric used as a stepsize between consecutive latitude lines. Value of zero
#' or `NULL` disables latitude lines.
#' @param region_args list of arguments passed to the plotting function for the region.
#'
#' @return NULL (invisibly)
#' @examples
#' ll <- seq(-20, 20, by = 5)
#' co <- expand.grid(lon = ll, lat = ll)
#' X <- s2pp(co, marks = co$lat)
#' plot(unmarx(X))
#' plot(X)
#' file.remove(img)
#' @export
plot.s2pp <- function(x, ..., add = FALSE, region = !add, longrid = 30, latgrid = 30,
                      region_args = list()){
  verifyclass(x, "s2pp")
  co <- s2coords(x)
  dots <- list(...)
  if(inherits(region, "s2region")){
    show.region <- TRUE
    ok <- s2contains(co, region)
  } else{
    if(is.logical(region)){
      show.region <- region
      region <- s2region(x)
      ok <- rep(TRUE, npoints(x))
    } else{
      stop("Argument ", sQuote("region"), " misspecified.")
    }
  }
  co <- co[ok,]
  if(is.marked(x)){
    marx <- marks(x, drop = FALSE)
    marx <- marx[ok, 1, drop=TRUE]
    if(is.numeric(marx)){
      cols <- spatstat.options("image.colfun")(256)
      cm <- colourmap(cols, range = range(marx))
      dots$col <- cm(marx)
    }
  }
  if(!add){
    globeearth(gdata = NULL, runlen = NULL, eye = dots$eye, top = dots$top)
    drawlonlatgrid(longrid = longrid, latgrid = latgrid, eye = dots$eye, top = dots$top)
  }
  if(show.region && !inherits(region, "s2")){
    do.call(plot, args = append(list(x = region, add = TRUE), region_args))
  }
  do.call(globepoints, args = append(list(loc = co), dots))
  return(invisible(NULL))
}

# Draws a lon,lat grid on globe
drawlonlatgrid <- function(longrid, latgrid, eye = NULL, top = NULL, ..., do.plot = TRUE, flat = FALSE){
  if(!is.null(latgrid) && latgrid!=0){
    if(length(latgrid)==1){
      # Make grid
      latgrid <- seq(-90, 90, by = latgrid)
      # Center the grid
      latgrid <- latgrid + (90-max(latgrid))/2
    }
    if(!flat)
      latgrid <- expand.grid(lon = c(seq(-180, 180, by=1), NA), lat = latgrid)
    if(do.plot){
      if(flat){
        abline(h = latgrid, ...)
      } else{
        globe::globelines(latgrid, eye = eye, top = top, ...)
      }
    }
  }
  if(!is.null(longrid) && longrid!=0){
    if(length(longrid)==1){
      # Make grid
      longrid <- seq(-180, 180, by = longrid)
      # Center the grid
      longrid <- longrid + (180-max(longrid))/2
    }
    if(!flat)
      longrid <- expand.grid(lat = c(seq(-90, 90, by=1), NA), lon = longrid)[,2:1]
    if(do.plot){
      if(flat){
        abline(v = longrid, ...)
      } else{
        globe::globelines(longrid, eye = eye, top = top, ...)
      }
    }
  }
  return(invisible(list(longrid = longrid, latgrid = latgrid)))
}


#' Plot the outline of a spherical polygon
#'
#' Plot the outline of a spherical polygon of class `"s2polygon"` in base
#' graphics via the [globe package][globe::globe-package].
#'
#' @param x Polygon of class `"s2polygon"` to plot, a list of loops on the
#' sphere or something interpretable as such.
#' @param eps Maximum distance between consecutive plotted vertices in the
#' polygon. For vertices separated by more than `eps` interpolating points will
#' be inserted. Default value corresponds to no interpolation.
#' @param ... parameters passed to [globe::globelines()] to control view point
#' (`eye`), orientation (`top`), and grapical parameters such as line color
#' (`col`), line width (`lwd`) etc.
#' @param add Logical to add outline of polygon to existing plot.
#' @param coast Logical to add the outline of the worlds coast lines.
#' @param flat Logical to plot on a flat earth map (simply projecting lon,lat
#' values to a rectangular region (-180,180)x(-90,90) using [globe::flatearth()]).
#' @inheritParams plot.s2pp
#'
#' @return NULL (invisibly)
#' @aliases plot.s2looplist plot.s2cell plot.s2cellid
#' @export
plot.s2polygon <- function(x, eps = pi * s2radius(x), ..., add = FALSE,
                           coast = FALSE, longrid = 30, latgrid = 30,
                           flat = FALSE){
  # Catch s2cells/s2cellids
  if(inherits(x, "s2cellid")){
    x <- s2cell(x)
  }
  if(inherits(x, "s2cell")){
    x <- x$vertices
  }

  if(!add){
    # Plot args for flat/globe-earth
    args <- resolve.defaults(list(col = "gray", lwd = 1), list(...))
    if(!coast)
      args <- append(args, list(gdata = NULL, runlen = NULL))
    plotfun <- if(flat) globe::flatearth else globe::globeearth
    do.call(plotfun, args)
    drawlonlatgrid(longrid = longrid, latgrid = latgrid, eye = args$eye, top = args$top, col = args$col, flat = flat)
  }
  eps <- eps/s2radius(x)
  loops <- s2looplist(x)
  if(eps<pi){
    loops <- lapply(loops, s2interpolate, eps = eps, radius = 1)
  }
  if(flat){
    s <- s2segments(loops)
    segments(s[,"fromlon"], s[,"fromlat"], s[,"tolon"], s[,"tolat"], ...)
  } else{
    lapply(loops, globelines, ...)
  }
  return(invisible(NULL))
}
#' @export
plot.s2looplist <- plot.s2polygon
#' @export
plot.s2cell <- plot.s2polygon
#' @export
plot.s2cellid <- plot.s2polygon

#' Start PNG plotting device for spherical objects
#'
#' This command starts a PNG plotting device and sets appropriate parameters for
#' plotting spherical objects. The device must be closed by the user by calling
#' `grDevices::dev.off()` or similar when no more graphics will be added.
#'
#' @param filename Name of output file. A temporary file will be generated if it
#' none is given.
#' @inheritDotParams grDevices::png
#' @param init Logical to initialize plot with a blank plot with
#' `xlim = c(-180,180)` and `ylim = c(-90,90)`.
#'
#' @return Name of output file.
#' @export
#'
#' @examples
#' img <- s2png()
#' globe::flatearth()
#' dev.off()
#' file.remove(img)
#'
s2png <- function(filename, ..., init = FALSE){
  if(missing(filename))
    filename <- tempfile(fileext = ".png")
  stopifnot(is.character(filename) && substr(filename, nchar(filename)-3, nchar(filename)) == ".png")
  args <- resolve.defaults(filename=filename, ..., list(width=2048, height=1024, bg="#AAEEFF", antialias="default", pointsize = 6))
  do.call.matched(png, arglist = as.list(args))
  par(mar = c(0,0,0,0),    pin = c(4,2),    pty = "m",    xaxs = "i",
      xaxt = "n",          xpd = FALSE,    yaxs = "i",    yaxt = "n")
  if(init)
    plot(0, 0, type = "n", axes = FALSE, xpd = FALSE, xlim = c(-180,180), ylim = c(-90,90))
  return(filename)
}

# s2png <- function(x, coast = FALSE, imgbg = "darkblue", ..., longrid = 30, latgrid = 30){
#   # Setup jpeg/png image and plot empty image
#   img <- s2png_init()
#   on.exit(dev.off())
#
#   # Plot background coastlines or lon,lat grid.
#   bgcol <- "darkgray"
#   bglwd <- 2
#   if(coast){
#     s <- globe::flatearth(do.plot = FALSE)
#     segments(s[,1], s[,2], s[,3], s[,4], col = bgcol, lwd = bglwd)
#   } else{
#     if(!is.null(latgrid) && latgrid!=0){
#       if(length(latgrid)==1){
#         latgrid <- seq(-90, 90, by = latgrid)
#         latgrid <- latgrid + (90-max(latgrid))/2
#       }
#       abline(h = latgrid, col = bgcol, lwd = bglwd)
#     }
#     if(!is.null(longrid) && longrid!=0){
#       if(length(longrid)==1){
#         longrid <- seq(-180, 180, by = longrid)
#         longrid <- longrid + (180-max(longrid))/2
#       }
#       abline(v = longrid, col = bgcol, lwd = bglwd)
#     }
#   }
#
#   if(inherits(x, "s2region")){
#     W <- x
#     P <- data.frame(lon = numeric(0), lat = numeric(0))
#   } else if(inherits(x, "s2pp")){
#     W <- s2region(x)
#     P <- as.data.frame.s2pp(x, format = "lon,lat")
#   } else{
#     stop("Can only handle s2region and s2pp.")
#   }
#   if(inherits(W, "s2polygon")){
#     s <- s2segments(W)
#     segments(s[,"fromlon"], s[,"fromlat"], s[,"tolon"], s[,"tolat"], ...)
#   }
#   points(P$lon, P$lat, ...)
#   return(img)
# }

s2interpolate <- function(x, eps = pi, radius = 1){
  if(eps>=pi*radius){
    return(x)
  }
  single <- FALSE
  if(is.data.frame(x) | is.matrix(x)){
    x <- list(x)
    single <- TRUE
  }
  # Put coords in right format (at the moment first 3d then lng,lat)
  x <- lapply(x, make_s2coords) ## Make 3d coords
  x <- lapply(x, globe::ensurelonlat)
  x <- lapply(x, as.data.frame)
  interpolate_one <- function(loop){
    lin <- s2::s2_make_line(loop$lon, loop$lat)
    len <- s2::s2_length(lin, radius = radius)
    dist <- c(seq(0, len, by = eps), len)
    pts <- s2::s2_interpolate(lin, distance = dist, radius = radius)
    coo <- wk::wk_coords(wk::wkb(s2_as_binary(pts)))
    return(data.frame(lon = coo$x, lat = coo$y))
  }
  rslt <- lapply(x, interpolate_one)
  if(single){
    rslt <- rslt[[1]]
  }
  return(rslt)
}

s2segments <- function(x, eps = pi, radius = 1){
  if(inherits(x, "s2polygon")){
    x <- x$loops
    radius <- s2radius(x)
  }
  if(inherits(x, "s2")){
    x <- list()
  }
  eps <- eps/radius
  loops <- x
  arcs <- NULL
  for(i in seq_along(loops)){
    l <- globe::ensure3d(loops[[i]], single = FALSE)
    if(eps<pi){
      l <- s2interpolate(l, eps = eps, radius = 1)
    }
    nn <- nrow(l)
    from <- globe::ensurelonlat(l[-nn,])
    to <- globe::ensurelonlat(l[-1,])
    fromto <- cbind(fromlat = from$lat, fromlon = from$lon, tolat = to$lat, tolon = to$lon)
    arcs <- rbind(arcs, fromto)
  }
  return(arcs)
}

s2coastlines <- function(flat = FALSE, ...){
  if(flat){
    s <- globe::flatearth(do.plot = FALSE)
  } else{
    s <- globe::globeearth(do.plot = FALSE)
  }
  segments(s[,1], s[,2], s[,3], s[,4], ...)
}

#' Interactive plotting function
#'
#' Generic to plot spherical data on an interactive sphere.
#'
#' @param x spherical data object to be plotted interactively.
#' @param ... Additional parameters passed to methods.
#' @export
iplot <- function(x, ...) {
  UseMethod("iplot")
}

#' Plot a spherical polygon on an interactive sphere With rotation and zoom
#'
#' Plot a spherical polygon on an interactive sphere using [threejs::globejs()].
#'
#' @param x Spherical polygon of class `"s2polygon"`, a list of loops on the
#' sphere or something interpretable as such.
#' @param ... Additional parameters passed to [threejs::globejs()] to control
#' the appearance of the plot. For instance `img` can be a file path to an image
#' to wrap on the sphere, and
#' `arcsLwd`, and `arcsColor` control the appearance of the plotted lines. When
#' `image = FALSE` the argument `arcsHeight` controls the height above the sphere
#' of the lines.
#' @param use_png Logical to plot the polygon on an png image which is then
#' wrapped on the sphere.
#' @param eps Maximum distance between consecutive plotted vertices in the
#' polygon. For vertices separated by more than `eps` interpolating points will
#' be inserted. Default value corresponds to no interpolation.
#'
#' @aliases iplot.s2looplist iplot.s2cell iplot.s2cellid
#' @export
iplot.s2polygon <- function(x, ..., use_png = FALSE, eps = pi * s2radius(x)){
  if(!requireNamespace("threejs"))
    stop("Package ", sQuote("threejs"), " is required for interactive plotting.")

  # Catch s2cells/s2cellids
  if(inherits(x, "s2cellid")){
    x <- s2cell(x)
  }
  if(inherits(x, "s2cell")){
    x <- x$vertices
  }
  args <- resolve.defaults(list(...), list(arcsHeight = 0.1, arcsLwd = 3, arcsColor = "red",
                           arcsOpacity = 1))
  if(!use_png){
    arcs <- s2segments(x, eps = eps)
    extraargs <- if(!is.null(arcs)) append(list(arcs = arcs), args) else args
    # img <- system.file("images/world.jpg", package = "threejs")
    globeplot <- do.call(threejs::globejs, extraargs)
    globeplot
  } else{
    givenimg <- args$img
    img <- s2png(init = !is.null(givenimg))
    if(!is.null(givenimg)){
      ext <- tail(strsplit(givenimg, ".", fixed = TRUE)[[1]], 1)
      switch(ext,
             jpg =,
             jpeg = {
               if(!requireNamespace("jpeg"))
                 stop("Need jpeg package to read jpeg image.")
               mat <- jpeg::readJPEG(givenimg)
               },
             png = {
               if(!requireNamespace("png"))
                 stop("Need png package to read png image.")
               mat <- png::readPNG(givenimg)
             },
             stop("Only png and jpeg/jpg extensions are recognised.")
      )
      rasterImage(mat, -180, -90, 180, 90)
    }
    plot(x, col = args$arcsColor, lwd = args$arcsLwd, flat = TRUE, add = !is.null(givenimg))
    dev.off()
    args <- resolve.defaults(list(img = img), ...)
    globeplot <- do.call(threejs::globejs, args)
    globeplot
  }
}

#' @export
iplot.s2looplist <- iplot.s2polygon

#' @export
iplot.s2cell <- iplot.s2polygon

#' @export
iplot.s2cellid <- iplot.s2polygon

#' Plot a spherical point pattern on an interactive sphere with rotation and zoom
#'
#' Plot a spherical point pattern on an interactive sphere using [threejs::globejs()].
#'
#' @param x Spherical point pattern of class `"s2pp"`.
#' @param ... Additional parameters passed to [iplot.s2polygon] and [threejs::globejs()] to control
#' the appearance of the plot. For instance `img` can be a file path to an image
#' to wrap on the sphere, and `color`, and `pointsize` control the appearance of
#' the plotted points.
#'
#' @export
iplot.s2pp <- function(x, ...){
  if(!requireNamespace("threejs"))
    stop("Package ", sQuote("threejs"), " is required for interactive plotting.")

  ll <- as.data.frame(x, format = "lon,lat")
  region <- s2region(x)
  if(inherits(region, "s2cap"))
    stop("Can't handle cap yet.")
  args <- list(lat = ll$lat, long = ll$lon, x = region)
  args <- resolve.defaults(args, list(...), list(color = "cyan", pointsize = 2, value = 0))
  globeplot <- do.call(iplot.s2polygon, args)
  globeplot
}
