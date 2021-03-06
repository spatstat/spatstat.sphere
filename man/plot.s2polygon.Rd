% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s2plotting.R
\name{plot.s2polygon}
\alias{plot.s2polygon}
\alias{plot.s2looplist}
\alias{plot.s2cell}
\alias{plot.s2cellid}
\title{Plot the outline of a spherical polygon}
\usage{
\method{plot}{s2polygon}(
  x,
  eps = pi * s2radius(x),
  ...,
  add = FALSE,
  coast = FALSE,
  longrid = 30,
  latgrid = 30,
  flat = FALSE
)
}
\arguments{
\item{x}{Polygon of class \code{"s2polygon"} to plot, a list of loops on the
sphere or something interpretable as such.}

\item{eps}{Maximum distance between consecutive plotted vertices in the
polygon. For vertices separated by more than \code{eps} interpolating points will
be inserted. Default value corresponds to no interpolation.}

\item{...}{parameters passed to \code{\link[globe:globelines]{globe::globelines()}} to control view point
(\code{eye}), orientation (\code{top}), and grapical parameters such as line color
(\code{col}), line width (\code{lwd}) etc.}

\item{add}{Logical to add outline of polygon to existing plot.}

\item{coast}{Logical to add the outline of the worlds coast lines.}

\item{longrid}{Numeric defining a grid of longitude lines to be plotted if
\code{add = FALSE}. Either a vector of numbers between -180 and 180 or a single
numeric used as a stepsize between consecutive longitude lines. Value of zero
or \code{NULL} disables longitude lines.}

\item{latgrid}{Numeric defining a grid of latitude lines to be plotted if
\code{add = FALSE}. Either a vector of numbers between -90 and 90 or a single
numeric used as a stepsize between consecutive latitude lines. Value of zero
or \code{NULL} disables latitude lines.}

\item{flat}{Logical to plot on a flat earth map (simply projecting lon,lat
values to a rectangular region (-180,180)x(-90,90) using \code{\link[globe:flatearth]{globe::flatearth()}}).}
}
\value{
NULL (invisibly)
}
\description{
Plot the outline of a spherical polygon of class \code{"s2polygon"} in base
graphics via the \link[globe:globe-package]{globe package}.
}
