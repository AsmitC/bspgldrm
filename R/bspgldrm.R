#' TItle: Fits a finite-support Dir-GLM model
#'
#' @import stats
#' @import mvtnorm
#' @import gldrm
#'
#' @param formula An object of class "formula".
#' @param data An optional data frame containing the variables in the model.
#' @param link Link function. Defaults to log.
#' Can be a character string to be passed to the
#' \code{make.link} function in the \code{stats} package (e.g. "identity",
#' "logit", or "log").
#' Alternatively, \code{link} can be a list containing three functions named
#' \code{linkfun}, \code{linkinv}, and \code{mu.eta}. The first is the link
#' function. The second is the inverse link function. The third is the derivative
#' of the inverse link function. All three functions must be vectorized.
#' @param bspgldrmControl Optional control arguments.
#' Passed as an object of class "bspgldrmControl", which is constructed by the
#' \code{bspgldrm.control} function.
#' See \code{bspgldrm.control} documentation for details.
#' @param thetaControl Optional control arguments for the theta update procedure.
#' Passed as an object of class "thetaControl", which is constructed by the
#' \code{theta.control} function.
#' See \code{theta.control} documentation for details.
#'
#' @return An S3 object of class "bspgldrm". See details.
#'
#' @details The arguments \code{linkfun}, \code{linkinv}, and \code{mu.eta}
#' mirror the "link-glm" class. Objects of this class can be created with the
#' \code{stats::make.link} function.
#'
#'The "bspgldrm" class is a list of the following items.
#' \itemize{
#' \item \code{samples} A list containing the MCMC samples for \code{f0} and \code{beta}
#' \item \code{formula} Model formula.
#' \item \code{data} Model data frame.
#' \item \code{link} Link function. If a character string was passed to the
#' \code{link} argument, then this will be an object of class "link-glm".
#' Otherwise, it will be the list of three functions passed to the \code{link} argument.
#' }
#'
#' @export
bspgldrm <- function(formula, data=NULL, link="log",
                     bspgldrmControl=bspgldrm.control(), thetaControl=theta.control())
{
  ## 1. Model Initialization
  mf <- stats::model.frame(formula, data)
  X  <- stats::model.matrix(attr(mf, "terms"), mf)
  attributes(X)[c("assign", "contrasts")] <- NULL
  y  <- stats::model.response(mf, type = "numeric")

  ## 2. Link Extraction
  if (is.character(link)) {
    link <- stats::make.link(link)
  } else if (!is.list(link) ||
             !all(c("linkfun", "linkinv", "mu.eta") %in% names(link))) {
    stop("link must be a string or a list containing linkfun, linkinv, mu.eta")
  }

  if (is.null(mu0)) {
    mu0 <- mean(y)
  } else if (mu0<=min(spt) || mu0>=max(spt)) {
    stop(paste0("mu0 must lie within the range of observed values. Choose a different ",
                "value or set mu0=NULL to use the default value, mean(y)."))
  }

  ## 3. MCMC Initialization
  betaStart <- bspgldrmControl$betaStart
  f0Start   <- bspgldrmControl$f0Start
  if (is.null(betaStart) || is.null(f0Start)) {
    gfit <- gldrm(formula       = formula,
                  data          = data,
                  link          = link,
                  mu0           = mu0,
                  gldrmControl  = gldrm.control(),
                  thetaControl  = thetaControl)
    if (is.null(betaStart)) betaStart <- gfit$beta
    if (is.null(f0Start))   f0Start   <- gfit$f0
  }

  init <- list(beta = betaStart,
               f0   = f0Start)
  spt <- sort(unique(y)) # Observed support

  ## 4. Call MCMC helper
  fit <- bspgldrmFit(
    X            = X,
    y            = y,
    spt          = spt,
    init         = init,
    joint.update = bspgldrmControl$joint.update,
    mu0          = mu0,
    rho          = bspgldrmControl$rho,
    burnin       = bspgldrmControl$burnin,
    thin         = bspgldrmControl$thin,
    save         = bspgldrmControl$save,
    thetaControl = thetaControl
  )

  ## 5. Output
  out <- list(
    samples = fit$samples,
    formula = formula,
    data    = data.frame(mf),
    link    = link
  )
  class(out) <- "bspgldrm"
  out
}
