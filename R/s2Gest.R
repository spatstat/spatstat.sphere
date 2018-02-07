#' G-function for spherical point pattern
#'
#' @param X Spherical point pattern of class `s2pp`
#' @param \dots Ignored.
#' @param r Optional. Vector of values for the argument \eqn{r} at which
#'   \eqn{G(r)} should be evaluated. Users are advised \emph{not} to specify
#'   this argument; there is a sensible default. If necessary, specify
#'   \code{rmax}.
#' @param rmax Optional. Maximum desired value of the argument \eqn{r}.
#' @param breaks This argument is for internal use only.
#' @param correction Optional. A character vector containing any selection of
#'   the options \code{"none"} and \code{"border"} (and possibly others in the
#'   future).
#' @return An object of class \code{"fv"}, see \code{\link{fv.object}}, which
#'   can be plotted directly using \code{\link{plot.fv}}.
#'
#'   Essentially a data frame containing columns \item{r}{the vector of values
#'   of the argument \eqn{r} at which the function \eqn{G} has been estimated }
#'   \item{theo}{the theoretical value \eqn{G(r) = 2 \pi (1-cos(r))}{G(r) = 2 *
#'   pi * (1 - cos(r))} for a stationary Poisson process (where r is in
#'   radians)} together with columns named \code{"none"}, and/or
#'   \code{"border"}, according to the selected edge corrections. These columns
#'   contain estimates of the function \eqn{G(r)} obtained by the edge
#'   corrections named.
#'
#' @examples
#' ## Random poins on earth
#' X <- s2runif(200, region = s2earth())
#' ## Estimate of G function up to distance 2000 km
#' G_X <- s2Gest(X, rmax = 2000)
#'
#' @export
s2Gest <- function(X, r = NULL, rmax = NULL, breaks = NULL, correction = NULL){
  angular_arg <- FALSE; degrees <- FALSE # This could potentially be an argument of the function in the long run.
  ratio <- FALSE # This could potentially be an argument of the function in the long run.

  # Check whether r values should be degrees (radians will be used and then changed to degrees at the end)
  dom <- s2region(X)
  radius <- s2radius(dom)

  # Figure out r-values
  rmaxdefault <- rmax %orifnull% (pi*radius)
  breaks <- handle.r.b.args(r, breaks, dom, pixeps = rmaxdefault/128, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max
  zeroes <- numeric(length(r))

  ## choose correction(s)
  if(is.null(correction)) {
    correction <- if(inherits(dom, "s2")) "none" else c("rs", "km")
  } else correction <- pickoption("correction", correction,
                                  c(none="none",
                                    border="rs",
                                    rs="rs",
                                    KM="km",
                                    km="km",
                                    Kaplan="km",
                                    best="km"),
                                  multi=TRUE)

  # Intensity
  areaW <- area(dom)
  n <- npoints(X)
  lambda <- n/areaW

  ## initialise output object
  theo <- 1 - exp( -lambda * 2*pi*radius^2 * (1-cos(r/radius)) )
  Odf <- data.frame(r = r, theo = theo)
  desc <- c(paste("distance argument", ifelse(!angular_arg, "r", "theta")),
            "theoretical uncorrected %s")
  ylab <- if(!angular_arg) quote(G(r)) else quote(G(theta))
  OO <- fv(Odf,
           argu = "r",
           ylab = ylab,
           valu = "theo",
           fmla = . ~ r,
           alim = c(0, ifelse(!degrees, rmax, rmax*360/(2*pi))),
           labl = c(ifelse(!angular_arg, "r", "theta"),
                    ifelse(!angular_arg, "%s[pois](r)", "%s[pois](theta)")),
           desc = desc,
           fname = "G",
           yexp = ylab
  )

  nnd <- nndist(X)

  if ("none" %in% correction) {
    ## Uncorrected for data on entire sphere.
    if (n <= 1)
      edf <- zeroes
    else {
      hh <- hist(nnd[nnd <= rmax], breaks = breaks$val,
                 plot = FALSE)$counts
      edf <- cumsum(hh)/length(nnd)
    }
    labl <- ifelse(!angular_arg, "hat(%s)(r)", "hat(%s)(theta)")
    OO <- bind.fv(OO, data.frame(raw = edf),
                  labl,
                  "uncorrected estimate of %s",
                  "raw")
  }

  if(any(correction %in% c("rs", "km"))) {
    ## calculate Kaplan-Meier and border correction (Reduced Sample) estimators
    if(n == 0)
      result <- data.frame(rs=zeroes, km=zeroes, hazard=zeroes, theohaz=zeroes)
    else {
      ##  distance to boundary
      bdry <- s2borderdist(X)
      ##  observations
      o <- pmin.int(nnd,bdry)
      ##  censoring indicators
      d <- (nnd <= bdry)
      result <- km.rs(o, bdry, d, breaks)
      result$theohaz <- 2 * pi * lambda * radius^2 * sin(r/radius)
      result <- as.data.frame(result[c("rs", "km", "hazard", "theohaz")])
    }
    ## add to fv object
    OO <- bind.fv(OO, result,
                  c("hat(%s)[bord](r)", "hat(%s)[km](r)",
                    "hat(h)[km](r)", "h[pois](r)"),
                  c("border corrected estimate of %s",
                    "Kaplan-Meier estimate of %s",
                    "Kaplan-Meier estimate of hazard function h(r)",
                    "theoretical Poisson hazard function h(r)"),
                  "km")

    ## modify recommended plot range
    attr(OO, "alim") <- range(r[result$km <= 0.9])
  }

  nama <- names(OO)
  fvnames(OO, ".") <- rev(setdiff(nama, c("r", "hazard", "theohaz")))
  unitname(OO) <- unitname(X)

  if(degrees) OO$r <- OO$r*360/(2*pi)

  return(OO)
}
