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
bspgldrm.control <- function(burnin=1000, thin=1, save=5000, rho=0.1, mu0=NULL,
                             betaStart=NULL, f0Start=NULL, joint.update=TRUE)
{
  if (burnin < 0)            stop("Number of burn-in samples must be >= 0")
  if (thin < 1)              stop("Thin must be >= 1")
  if (save < 1)              stop("Must save at least one iteration")
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

#' Main MCMC function
#'
#' This function is called by the main \code{bspgldrm} function.
#'
#' @keywords internal
bspgldrmFit <- function(X, y,
                        linkfun, linkinv, mu.eta,
                        control, thetaControl)
{
  ## 1. Extract control parameters
  burnin       <- control$burnin
  thin         <- control$thin
  save         <- control$save
  rho          <- control$rho
  mu0          <- control$mu0
  joint.update <- control$joint.update
  betaStart    <- control$betaStart
  f0Start      <- control$f0Start

  spt <- sort(unique(y)) ## Observed support

  ## 2. Initialize mu0 if not provided by the user
  if (is.null(mu0)) {
    mu0 <- mean(y)
  } else if (mu0 <= min(spt) || mu0 >= max(spt)) {
    stop(paste0("mu0 must lie within the range of observed values. Choose a different ",
                "value or set mu0=NULL to use the default value, mean(y)."))
  }

  ## 3. MCMC Initialization
  if (is.null(betaStart) || is.null(f0Start)) {
    ### Fit gldrm by default
    gfit <- gldrm(formula      = formula,
                  data         = data,
                  link         = link,
                  mu0          = mu0,
                  gldrmControl = gldrm.control(),
                  thetaControl = thetaControl)
    if (is.null(betaStart)) betaStart <- gfit$beta
    if (is.null(f0Start))   f0Start   <- gfit$f0
  }

  init <- list(beta = betaStart, f0 = f0Start)
  X    <- as.matrix(X)
  n    <- length(y)
  p    <- ncol(X)
  l    <- length(spt)
  iter <- burnin + thin * save

  ### 3.1 Beta and f0
  beta_samples <- matrix(NA, nrow = save, ncol = p)
  f0_samples   <- matrix(NA, nrow = save, ncol = l)

  beta <- init$beta
  f0   <- init$f0
  beta_samples[1, ] <- beta
  f0_samples[1, ]   <- f0

  ### 3.2 Dirichlet prior
  ind_mt      <- outer(y, spt, `==`) * 1
  alpha       <- 1
  dir_pr_parm <- alpha * colMeans(ind_mt)
  eps         <- 1e-6
  dir_pr_parm <- dir_pr_parm + eps

  # 3.3 Theta
  mu      <- linkinv(X %*% beta)                    # Updated for general link
  out     <- tht_sol(spt, f0, mu, NULL)
  tht     <- out$tht
  btht    <- out$btht
  bpr2    <- out$bpr2
  f0_y    <- f0y(y, spt, f0)

  ## 4. MCMC loop
  for (r in 2:iter) {
    ### 4.1 Beta update
    Sig <- Sigma_beta(X, mu, bpr2, rho, linkfun, mu.eta)
    b_out <- ifelse(joint.update,
                    beta_update_joint(X, y, spt, beta, Sig, f0,
                                      tht, bpr2, btht, rho, linkinv),
                    beta_update_separate(X, y, spt, beta, Sig, f0,
                                         tht, bpr2, btht, rho, linkinv))
    beta <- b_out$cr_bt
    tht  <- b_out$cr_tht
    btht <- b_out$cr_btht
    bpr2 <- b_out$cr_bpr2
    mu   <- linkinv(X %*% beta)                   # Updated for general link

    ### 4.2 f0 update
    propsl_dir_parm <- dir_parm(y, tht, btht, dir_pr_parm, ind_mt)
    out             <- f0_update(y, spt, f0, f0_y, propsl_dir_parm,
                                 mu, tht, bpr2, btht, dir_pr_parm, ind_mt)
    f0    <- out$cr_f0
    f0_y  <- out$cr_f0y
    tht   <- out$cr_tht
    btht  <- out$cr_btht
    bpr2  <- out$cr_bpr2

    # 4.3 Storage
    if (r > burnin & r %% thin == 0) {
      j <- (r - burnin) / thin
      beta_samples[j, ] <- beta
      f0_samples[j, ]   <- f0
    }
  }

  ## 5. Output
  list(samples = list(beta = beta_samples,
                      f0   = f0_samples))
}

