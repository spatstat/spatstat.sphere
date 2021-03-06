% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/s2region.R
\name{s2lonlatbox}
\alias{s2lonlatbox}
\alias{s2lonlatbox.s2region}
\alias{s2lonlatbox.default}
\title{Generic function to extract or construct a box in longitude and latitude coordinates}
\usage{
s2lonlatbox(...)

\method{s2lonlatbox}{s2region}(x, ...)

\method{s2lonlatbox}{default}(lon = c(-180, 180), lat = c(-90, 90), ...)
}
\arguments{
\item{...}{Parameters passed to methods. For \code{s2lonlatbox.s2region} only
first argument (named \code{x}) is used and the rest ignored. For
\code{s2lonlatbox.default} the first two arguments (named \code{lon} and \code{lat}) are
used to define the region and the rest are passed to \code{\link[=s2]{s2()}} to define the
radius and units of the underlying sphere.}

\item{x}{Object of class \code{"s2region"}.}

\item{lon}{Numeric of length 2 specifying longitude interval in degrees. For
boxes not crossing the 180 degree meridian the input should be such that
\verb{-180 <= lon[1] < lon[2] <= 180}. For boxes crossing the 180 degree meridian
the input shoud be such that \verb{-180 <= lon[2] < lon[1] <= 180}.}

\item{lat}{Numeric of length 2 with \verb{-90 <= lat[1] < lat[2] <= 90} specifying
latitude interval in degrees.}
}
\value{
Box in longitude and latitude coordinates (of class \code{"s2lonlatbox"}).
}
\description{
Generic function to extract or construct a box in longitude and latitude coordinates.
}
\section{Methods (by class)}{
\itemize{
\item \code{s2region}: Extraction method for \code{"s2region"}.

\item \code{default}: Default method to generate \code{"s2lonlatbox"} from a
longitude and latitude interval.
}}

\examples{
cap <- s2cap(c(1, 0, 0), h = 0.1)
s2lonlatbox(cap)
s2lonlatbox(lon = c(-30, 30), lat = c(45, 90))
## Very large box crossing the 180th meridian
big <- s2lonlatbox(lon = c(1, -1), lat = c(-90, 90))
area(big)
## Very small box **not** crossing the 180th meridian
small <- s2lonlatbox(lon = c(-1, 1), lat = c(-90, 90))
area(small)

}
\seealso{
s2lonlatbox.s2region s2lonlatbox.default
}
