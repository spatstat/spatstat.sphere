#' K-function for spherical point pattern
#'
#' @param X Spherical point pattern of class `s2pp`
#' @param \dots Ignored.
#' @param r Optional. Vector of values for the argument \eqn{r} at which
#'   \eqn{K(r)} should be evaluated. Users are advised \emph{not} to specify
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
#'   of the argument \eqn{r} at which the function \eqn{K} has been estimated }
#'   \item{theo}{the theoretical value \eqn{K(r) = 2 \pi (1-cos(r))}{K(r) = 2 *
#'   pi * (1 - cos(r))} for a stationary Poisson process (where r is in
#'   radians)} together with columns named \code{"none"}, and/or
#'   \code{"border"}, according to the selected edge corrections. These columns
#'   contain estimates of the function \eqn{K(r)} obtained by the edge
#'   corrections named.
#' @export
s2Kest <- function(X, ..., r = NULL, rmax = NULL, breaks = NULL, correction = NULL){
  angular_arg <- FALSE; degrees <- FALSE # This could potentially be an argument of the function in the long run.
  ratio <- FALSE # This could potentially be an argument of the function in the long run.

  # Check whether r values should be degrees (radians will be used and the chaged to degrees at the end)
  dom <- s2region(X)
  radius <- s2radius(dom)

  # Figure out r-values
  rmaxdefault <- rmax %orifnull% (pi*radius)
  breaks <- handle.r.b.args(r, breaks, dom, pixeps = rmaxdefault/128, rmaxdefault=rmaxdefault)
  r <- breaks$r
  rmax <- breaks$max

  # Choose correction based on domain if unset by user
  best_correction <- ifelse(inherits(dom, "s2"), "none", "border")
  if(missing(correction) || is.null(correction)){
    correction <- best_correction
  }
  correction <- pickoption("correction", correction,
                           c(none="none", border="border", best = best_correction),
                           multi=TRUE)

  # Intensity
  areaW <- area(dom)
  n <- npoints(X)
  lambda <- n/areaW
  lambda2 <- lambda * (n-1)/areaW

  ## initialise output object
  theo <- 2*pi*radius^2*(1-cos(r/radius))
  Odf <- data.frame(r = r, theo = theo)
  desc <- c(paste("distance argument", ifelse(!angular_arg, "r", "theta")),
            "theoretical uncorrected %s")
  ylab <- if(!angular_arg) quote(K(r)) else quote(K(theta))
  OO <- fv(Odf,
           argu = "r",
           ylab = ylab,
           valu = "theo",
           fmla = . ~ r,
           alim = c(0, ifelse(!degrees, rmax, rmax*360/(2*pi))),
           labl = c(ifelse(!angular_arg, "r", "theta"),
                    ifelse(!angular_arg, "%s[pois](r)", "%s[pois](theta)")),
           desc = desc,
           fname = "K",
           unitname = s2radius(dom, format = "list")$units,
           yexp = ylab
  )

  close <- closepairs(X, rmax = rmax)
  delta <- close$d

  if(any(correction == "none")) {
    ## Uncorrected for data on entire sphere.
    hh <- hist(delta, breaks = c(r, Inf), plot = FALSE)$counts
    khat <- c(0, cumsum(hh[-length(hh)]))
    khat <- khat/(lambda2*areaW)

    ## uncorrected estimate
    labl <- ifelse(!angular_arg, "hat(%s)(r)", "hat(%s)(theta)")
    OO <- bind.fv(OO, data.frame(un=khat),
                  labl,
                  "uncorrected estimate of %s",
                  "un")
  }

  if(any(correction == "border")) {
    ## border method
    ## Compute distances to boundary
    b <- s2borderdist(X)
    I <- close$i
    bI <- b[I]
    ## apply reduced sample algorithm
    RS <- Kount(close$d, bI, b, breaks)
    numKb <- RS$numerator
    denKb <- lambda * RS$denom.count
    labl <- ifelse(!angular_arg, "hat(%s)[bord](r)", "hat(%s)[bord](theta)")
    OO <- bind.ratfv(OO,
                     data.frame(border=numKb),
                     data.frame(border=denKb),
                     labl,
                     "border-corrected estimate of %s",
                     "border",
                     ratio=ratio)
  }

  if(degrees) OO$r <- OO$r*360/(2*pi)

  return(OO)
}
