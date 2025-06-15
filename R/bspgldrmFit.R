#' Control arguments for the \code{bspgldrm} algorithm.
#'
#' This function returns control arguments for the \code{bspgldrm} algorithm.
#' Each argument has a default value, which will be used unless a different
#' value is provided by the user.
#'
#' @param burnin Number of burn-in MCMC iterations. Defaults to 1000.
#' @param thin Factor by which to thin MCMC iterations. Defaults to 1.
#' @param save Number of MCMC samples to return. Defaults to 5000.
#' @param rho MCMC update step size, a (scalar) in (0, 1]. Defaults to 1.
#' @param mu0 Mean of the reference distribution. The reference distribution is
#' not unique unless its mean is restricted to a specific value. This value can
#' be any number within the range of observed values, but values near the boundary
#' may cause numerical instability. This is an optional argument with \code{mean(y)}
#' being the default value.
#' @param betaStart Initial value for the regression coefficients \code{beta}.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param f0Start Initial value for the reference distribution \code{f0}.
#' Defaults to the output obtained by fitting \code{gldrm}.
#' @param joint.update Logical specifying whether to update beta jointly or separately.
#' Defaults to \code{TRUE}
#'
#' @return Object of S3 class "bspgldrmControl"
#'
#' @export
bspgldrm.control <- function(burnin=1000, thin=1, save=5000, rho=1, mu0=NULL,
                             betaStart=NULL, f0Start=NULL, joint.update=TRUE)
{

  if (burnin < 0) stop("Number of burn-in samples must be >= 0")
  if (thin < 1) stop("Thin must be >= 1")
  if (save < 1) stop("Must save at least one iteration")
  if (!(rho <= 1 & rho > 0)) stop("rho must lie in (0, 1]")
  ctrl <- list(burnin      = burnin,
               thin        = thin,
               save        = save,
               rho.        = rho,
               mu0         = mu0,
               betaStart   = betaStart,
               f0Start     = f0Start,
               joint.update= joint.update,)
  class(ctrl) <- "bspgldrmControl"
  ctrl
}


