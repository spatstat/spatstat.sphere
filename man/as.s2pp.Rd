% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s2pp.R
\name{as.s2pp}
\alias{as.s2pp}
\alias{as.s2pp.s2pp}
\alias{as.s2pp.data.frame}
\alias{as.s2pp.matrix}
\alias{as.s2pp.default}
\title{Convert data to class s2pp}
\usage{
as.s2pp(X, ...)

\method{as.s2pp}{s2pp}(X, ...)

\method{as.s2pp}{data.frame}(X, ..., region = NULL)

\method{as.s2pp}{matrix}(X, ..., region = NULL)

\method{as.s2pp}{default}(X, ..., region = NULL)
}
\arguments{
\item{X}{Data which will be converted into a point pattern}

\item{...}{
  Arguments passed on to \code{\link[=s2]{s2}}
  \describe{
    \item{\code{radius}}{Numeric. Defaults to 1. For convenience this can also be an
object of class \code{"s2region"} from where the radius and unitname then will be
extracted.}
    \item{\code{unitname}}{Optional. Name of unit of length. Either a single character
string, or a vector of two character strings giving the singular and plural
forms, respectively.}
  }}

\item{region}{Region on the sphere where the points occur. Object of class
\code{"s2region"}.}
}
\value{
An object of class \code{"s2pp"} (see \code{\link{s2pp}}) describing
the spherical point pattern and its region.
}
\description{
Tries to coerce any reasonable kind of data to a spherical point pattern (an
object of class \code{"s2pp"}) for use by the \pkg{spatstat.sphere} package).
}
\details{
Converts the dataset \code{X} to a point pattern (an object of class
\code{"s2pp"}).

The function \code{as.s2pp} is generic, with methods for the classes
\code{"s2pp"}, \code{"matrix"}, \code{"data.frame"} and a default method.

Point pattern datasets can also be created by the function \code{\link{s2pp}}.
}
\section{Methods (by class)}{
\itemize{
\item \code{as.s2pp(s2pp)}: Method for s2pp objects.

\item \code{as.s2pp(data.frame)}: Method for data frames.

\item \code{as.s2pp(matrix)}: Method for matrices.

\item \code{as.s2pp(default)}: Default method.

}}
\keyword{manip}
\keyword{spatial}
