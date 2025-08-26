% File c:/sae2/man/mvrnormSeries.Rd
\name{mvrnormSeries}
\alias{mvrnormSeries}
\title{Generate data under the Dynamic or Rao-Yu Time Series Models}
\description{
  Function to generate data under a Rao-Yu time series model, a
  dynamic model, or a mixture of both. The function can produce either
  univariate or multivariate observations. All components of the returned
  random variable have unconditional mean zero. The function calls 
  \code{mvrnorm} in \pkg{MASS}.
}
\usage{
mvrnormSeries(NV=1,  D,  TI,  sigma.e, 
   rho.dyn,  sigma.v.dyn,  sigma.u.dyn, rho.u.dyn, 
   rho.RY,   sigma.v.RY,   sigma.u.RY,   rho.u.RY,  tol=1e-6,
   conditional.mean=FALSE)
}
\arguments{
  \item{NV}{The number of variables.}
  \item{D}{The number of domains.}
  \item{TI}{The number of time instances (constant for all domains).}
  \item{sigma.e}{The covariance matrix for the variation due to sampling,
    specified either as a single square matrix with \code{NV*TI} rows and 
    columns, or as a list of \code{D} matrices, each with \code{NV*TI} rows 
    and columns. The rows should vary over \code{TI} more quickly than over
    \code{NV}. Sampling covariances between domains are assumed to be zero.
    }
  \item{rho.dyn}{The temporal correlation parameter in the dynamic model.
    This parameter is not a true correlation, however, and it may exceed 1.}
  \item{sigma.v.dyn}{A vector of length \code{NV} with the v (nu) component
    of the variance under the dynamic model. This component measures the
    variance of the random effect at the beginning of the series.}
  \item{sigma.u.dyn}{A vector of length \code{NV} with the u component
    of the variance under the dynamic model.}
  \item{rho.u.dyn}{For \code{NV>1}, the cross-sectional correlation(s) under
    the model. The vector should specify \code{(NV*(NV-1))/2} correlations, in 
    lower triangular form. For example, for \code{NV=3}, the 
    correlations should be specified in the order (2,1), (3,1), (3,2).}
  \item{rho.RY}{The temporal correlation parameter in the Rao-Yu model. This
    is a true correlation, unlike the parameter in the dynamic model.}
  \item{sigma.v.RY}{A vector of length \code{NV} with the v (nu) component
    of the variance under the Rao-Yu model. This component reflects a constant
    random effect for each domain unchanging over time.}
  \item{sigma.u.RY}{A vector of length \code{NV} with the u component
    of the variance under the Rao-Yu model.}
  \item{rho.u.RY}{For \code{NV>1}, the cross-sectional correlation under
    the model. The vector should specify \code{(NV*(NV-1))/2} correlations, in 
    lower triangular form. For example, for \code{NV=3}, the 
    correlations should be specified in the order (2,1), (3,1), (3,2).}
  \item{tol}{A tolerance parameter used by \code{mvrnorm} in \pkg{MASS} to
    determine if the covariance matrix is non-singular.}
  \item{conditional.mean}{If true, the function will also return the 
    generated values of the random effects.}
}
\details{
The function assembles the covariance matrix from the covariance matrix 
under the dynamic model (if specified), the Rao-Yu model (if specified)
and a required sampling covariance matrix.
}
\value{If \code{conditional.mean=FALSE}, then for \code{NV=1}, a multivariate 
       normal random vector with mean zero and length \code{D*TI}. 
       For \code{NV>1}, a matrix with \code{D*TI} rows and \code{NV} 
       columns.
       
       If \code{conditional.mean=TRUE}, a list with the first element as above 
       and a second element that is the sum of the random effects without
       the sampling error. Simulation studies can evaluate the small area
       estimates using the first element of the list as input against
       the second element of the list, which is the target of the small area
       estimation.            
}
\author{
Robert E. Fay
}
\examples{
set.seed(7)
mvrnormSeries(D=2, TI=5, sigma.e=diag(5), rho.dyn=.8,
    sigma.v.dyn=2, sigma.u.dyn=.72, conditional.mean=TRUE)
mvrnormSeries(NV=2, D=2, TI=5, sigma.e=diag(10), rho.dyn=.8,
    sigma.v.dyn=2, sigma.u.dyn=.72, rho.u.dyn=.8)
}
