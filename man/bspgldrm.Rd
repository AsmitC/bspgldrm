% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bspgldrm.R
\name{bspgldrm}
\alias{bspgldrm}
\title{Fits a finite-support Dir-GLM model}
\usage{
bspgldrm(
  formula,
  data = NULL,
  link = "log",
  mb = NULL,
  sb = NULL,
  dir_pr_parm = NULL,
  bspgldrmControl = bspgldrm.control(),
  thetaControl = theta.control()
)
}
\arguments{
\item{formula}{An object of class "formula".}

\item{data}{An optional data frame containing the variables in the model.}

\item{link}{Link function. Defaults to log.
Can be a character string to be passed to the
\code{make.link} function in the \code{stats} package (e.g. "identity",
"logit", or "log").
Alternatively, \code{link} can be a list containing three functions named
\code{linkfun}, \code{linkinv}, and \code{mu.eta}. The first is the link
function. The second is the inverse link function. The third is the derivative
of the inverse link function. All three functions must be vectorized.}

\item{mb}{Prior mean for beta. Defaults to a p-length vector whose entries are all 0.}

\item{sb}{Vector containing the diagonal entries of the prior variance-covariance
matrix for beta. Defaults to the p-dimensional identity matrix.}

\item{dir_pr_parm}{Dirichlet prior parameter for f0. Defaults to the observed
response frequency distribution. If specified, it should be a p-length vector
with positive entries.}

\item{bspgldrmControl}{Optional control arguments.
Passed as an object of class "bspgldrmControl", which is constructed by the
\code{bspgldrm.control} function.
See \code{bspgldrm.control} documentation for details.}

\item{thetaControl}{Optional control arguments for the theta update procedure.
Passed as an object of class "thetaControl", which is constructed by the
\code{theta.control} function.
See \code{theta.control} documentation for details.}
}
\value{
An S3 object of class "bspgldrm". See details.
}
\description{
Fits a finite-support Dir-GLM model
}
\details{
The arguments \code{linkfun}, \code{linkinv}, and \code{mu.eta}
mirror the "link-glm" class. Objects of this class can be created with the
\code{stats::make.link} function.

The "bspgldrm" class is a list of the following items.
\itemize{
\item \code{samples} A list containing the MCMC samples for \code{f0} and \code{beta}.
\item \code{mb} Prior mean for beta.
\item \code{sb} Diagonal entries of the prior variance-covariance matrix for beta.
\item \code{dir_pr_parm} Dirichlet prior parameter.
\item \code{formula} Model formula.
\item \code{data} Model data frame.
\item \code{link} Link function. If a character string was passed to the
\code{link} argument, then this will be an object of class "link-glm".
Otherwise, it will be the list of three functions passed to the \code{link} argument.
}
}
\examples{
data(iris, package="datasets")

# Fit a bspgldrm with log link
fit <- bspgldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
                data=iris)
fit

# Fit a bspgldrm with custom link function
link <- list()
link$linkfun <- function(mu) log(mu)^3
link$linkinv <- function(eta) exp(eta^(1/3))
link$mu.eta <- function(eta) exp(eta^(1/3)) * 1/3 * eta^(-2/3)
fit2 <- bspgldrm(Sepal.Length ~ Sepal.Width + Petal.Length + Petal.Width + Species,
                 data=iris, link=link)
fit2

}
