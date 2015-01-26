% File ~/sae2/man/remldyn.R
\name{reml.dyn}
\alias{reml.dyn}
\alias{mle.dyn}
\alias{reml.Rao.Yu}
\alias{mle.Rao.Yu}
\title{Internal fitting functions for Dynamic and Rao-Yu models}
\description{
  Functions called by \code{eblupDyn} or \code{eblupRY} to produce EBLUP 
  small area estimates of the dynamic or Rao-Yu time series models through
  either ML or REML estimation of the variance components. The functions 
  can fit univariate or multivariate models.
}
\usage{
reml.dyn(y, X, M, T, NV=1, vcov_e, maxiter=100, iter.tol=.1e-5,
         sig2_u = 1, sig2_v=1, rho=.8, rho_u =.4, delta=NULL, 
         rho.fixed=NULL, y.include=NULL, ids=NULL, contrast.matrix=NULL,  
         baby.steps=TRUE, dampening=1, iter.history=FALSE,
         sig2.min.factor=.0001, max.rho_u=.98, tol=.Machine$double.eps,
         y.rescale=NULL) 

mle.dyn(y, X, M, T, NV=1, vcov_e, maxiter=100, iter.tol=.1e-5, 
         sig2_u = 1, sig2_v=1, rho=.8, rho_u =.4, delta=NULL,
         rho.fixed=NULL, y.include=NULL, ids=NULL, contrast.matrix=NULL,
         baby.steps=TRUE, dampening=1, iter.history=FALSE,
         sig2.min.factor=.0001, max.rho_u=.98, tol=.Machine$double.eps,
         y.rescale=NULL) 

reml.Rao.Yu(y, X, M, T, NV=1, vcov_e, maxiter=100, iter.tol=.1e-5, 
         sig2_u = 1, sig2_v=1, rho=.8, rho_u =.4, delta=NULL,
         rho.fixed=NULL, y.include=NULL, ids=NULL, contrast.matrix=NULL,    
         baby.steps=TRUE, dampening=.9, iter.history=FALSE,
         sig2.min.factor=.0001, max.rho_u=.98, max.rho=.98,
         tol=.Machine$double.eps, y.rescale=NULL)

mle.Rao.Yu(y, X, M, T, NV=1, vcov_e, maxiter=100, iter.tol=.1e-5,  
         sig2_u = 1, sig2_v=1, rho=.8, rho_u =.4, delta=NULL, 
         rho.fixed=NULL, y.include=NULL, ids=NULL, contrast.matrix=NULL, 
         baby.steps=TRUE, dampening=.9, iter.history=FALSE,
         sig2.min.factor=.0001, max.rho_u=.98, max.rho=.98,
         tol=.Machine$double.eps, y.rescale=NULL)

}
\arguments{
  \item{y}{For a univariate model, the dependent variable sorted in
    ascending order by time within domain. For a multivariate model, the 
    dependent variables sorted in ascending order by time within variable
    within domain.}
  \item{X}{A matrix of independent variables with the same number of 
    rows as the length of \code{y}.}
  \item{M}{The total number of domains, equivalent to \code{D} in 
      \code{eblupDyn} and \code{eblupRY}.}
  \item{T}{The number of time instances (constant for all domains).}
  \item{NV}{The number of dependent variables.}
  \item{vcov_e}{For the univariate model, the sampling covariance matrix 
    for the direct estimates of the \code{M*T} elements of the dependent 
    variable. The covariance matrix should be in the form of a square
    matrix with \code{M*T} rows and columns. Non-zero covariances between 
    domains are not allowed, so the matrix must have block
    diagonal form with \code{M} blocks, each of which is a square matrix
    with \code{T} rows and columns. Note that within domain, non-zero 
    covariances are allowed over time. 
    
    For the multivariate model, the square covariance matrix for the 
    \code{M*NV*T} the dependent variables. The matrix should be in the
    form of a square matrix with \code{M*NV*T} rows and columns. Time 
    should vary within variable, which should vary within domain.
    Non-zero covariances between domains are not allowed, but non-zero 
    covariances may be present across time and between components.}
  \item{maxiter}{The maximum number of iterations allowed for the 
    Fisher-scoring algorithm, with a default value of 100.}
  \item{iter.tol}{The convergence tolerance limit for the 
    Fisher-scoring algorithm, with a default value of .000001.}
  \item{sig2_u}{An initial starting value or values for the variance of the
    random increments.}
  \item{sig2_v}{An initial starting value or values for a domain level
    random effect. In the Rao-Yu model, the random effect is constant over
    time, whereas in the Dynamic model it is an initial effect 
    subject to dampening over time.}
  \item{rho}{The correlation across time. This correlation is assumed to
    be the same for the dependent variables in the multivariate model.}
  \item{rho_u}{For NV>1 only, the (NV*(NV-1))/2 starting values for
    the correlations between the random effects of the different dependent
    variables. If a single value is given, it will be used for 
    (NV*(NV+1))/2 components. The sort order corresponds to a lower triangle
    of the covariance matrix.}
  \item{delta}{The random effect components in the preferred internal order.
    Specification of \code{delta} will override any specification of \code{rho},
    \code{sig2_u}, \code{sig2_v}, or \code{rho_u}.}
  \item{rho.fixed}{If TRUE, the value of \code{rho} imbedded in \code{delta},
    if specified, or else given by \code{rho} will remain fixed during the
    iterations. Among other features, this allows the likelihood function
    for trial values of \code{rho} to be computed at the maximum over the other
    random effect parameters.}
  \item{y.include}{If specified, vector of length \code{M} to indicate 
    which domains to include in the estimation, with 1 signalling inclusion
    and 0 exclusion. Estimates for the excluded domain will be based on
    the fixed effect model only.}
  \item{ids}{A data frame with \code{M} rows giving ids for each of the
             domains. The data frame is copied to the returned object.}
  \item{contrast.matrix}{A matrix of coefficients of contrasts. The matrix
             must have \code{T*NV} rows, but it can contain an arbitrary number
             of columnts. Within each domain, the coefficients are applied
             to the \code{T*NV} EBLUP estimates.}
  \item{baby.steps}{Unless specified as FALSE, the first five iterations of the
             Fisher scoring algorithm are dampened by factors of 
             \code{c(.0625, .125, .25, .5, .75)}. These heuristically derived 
             factors appear to lessen drastic overshooting of the true maximum
             in the initial iterations.}
  \item{dampening}{A factor used to dampen the changes in the random
             effect parameters. Unlike \code{baby.steps}, its effect persists during
             all of the iterations until convergence. Note that the "factory
             setting" of this parameter is 1 for the dynamic model but .9 for 
             the Rao-Yu model.}
  \item{iter.history}{If TRUE, key values are saved during each of the 
             iterations and included as additional items in the returned list:
     \itemize{
       \item \code{delta.hist:} values of \code{delta} at each step.
       \item \code{likelihood.hist:} values of the log-likelihood (for ML) or
                restricted log-likelihood (for REML) at each step.
       \item \code{adj.hist:} adjustment at each step to prevent violation of
                boundary constraints.
       \item \code{s:} value of \code{s} used to compute next step in the
                iteration.
       \item \code{ix.hist:}  components of \code{delta} subject to alteration
                in this step.
       \item \code{adj.factor:} the final multiplier of the iteration step.
       \item \code{warning.hist:} an indicator of warning messages generated
                at each step.
       }
     }
  \item{sig2.min.factor}{A factor to multiply the minimum direct variance to
         set a minimum value for any of the variance components.}
  \item{max.rho_u}{A maximum allowed value for the estimated \code{rho_u}.}
  \item{tol}{A tolerance value used by matrix routines to prevent numerical
         instability. The value may be set to a lower value to encourage
         covergence, but appropriate caution should be applied.}
  \item{y.rescale}{A scaler multiplier for all of the \code{y} values. If the
         \code{y} values are either too small or two large, the information 
         matrix may become unstable. Setting this value to 1 has no effect; 
         setting it to 10 or 100 rescales very small \code{y} values to a more
         appropriate range. Similarly, positive values less than 1 may be used
         to rescale large \code{y} values. The effect of rescaling is removed 
         before normal return from the function, within the limits of normal 
         precision.}
  \item{max.rho}{For the Rao-Yu model only, a maximum allowed value for 
         \code{rho}.}
    
}
\details{
A typical model has the form response ~ terms where response is the (numeric) response vector and terms is a series of terms which specifies a linear predictor for response. A terms specification of the form first + second indicates all the terms in first together with all the terms in second with duplicates removed. 

A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} for more details of allowed formulae. 

The additional parameters include \code{contrast.matrix}, which specifies linear combinations of estimates within domains, such as the sum over components or across time. Corresponding MSE estimates are provided.

Other parameters affect convergence or provide starting values.
}
\value{
   \item{eblup}{In the univariate case, a vector of length \code{M*T} with 
   the eblup estimates. In the multivariate case, a matrix of M*T rows 
   and NV columns.}
   \item{fit}{A list summarizing the fit of the model with the following:
     \itemize{
       \item \code{model:} form of the model: T - Dynamic or RaoYu; REML
            or ML.
       \item \code{covergence:} a logical value indicating whether the 
                 convergence criterion was met.
      }
   }
   \item{parm}{A labelled vector with the estimated variance components, 
         correlations, and number of iterations.}
   \item{coef}{A labelled vector of coefficients of the model or models.}
   \item{ids}{A data frame with \code{D} rows and one or more columns of 
         numeric or character domain identifiers.}
   \item{delta}{An ordered vector of the variance components, which may be 
         used as starting values for additional iterations.}
   \item{eblup.mse}{MSE estimates for eblup.}
   \item{eblup.g1}{The g1 term of the MSE estimate.}
   \item{eblup.g2}{The g2 term of the MSE estimate.}
   \item{eblup.g3}{The g3 term of the MSE estimate.}
   \item{est.fixed}{Estimates based on fixed effects only.}
   \item{est.fixed.var}{The variance-covariance matrix for the estimates in 
         \code{coef}.}
   \item{eblup.wt1}{Weights given to the direct estimate in forming \code{eblup}.}
   \item{eblup.wt2}{Weights given to the direct estimate, including
        effects through estimating the fixed effect coefficients.}
   \item{contrast.est}{Estimates requested by the specified contrasts.}
   \item{contrast.mse}{MSE estimates for \code{contrast.est}.}
   \item{contrast.g1}{The g1 term in the estimation of \code{contrast.mse}.}
   \item{contrast.g2}{The g2 term in the estimation of \code{contrast.mse}.}
   \item{contrast.g3}{The g3 term in the estimation of \code{contrast.mse}.}
   \item{contrast.fixed.est}{Contrast estimates based on the fixed effect 
             model.}
   \item{contrast.fixed.var}{Variance estimates for the fixed effect model.}
   \item{contrast.wt1}{Weight wt1 given to the direct estimate in estimating 
       the contrasts.}
   \item{contrast.wt2}{Weight wt2 in estimating the contrasts.}
   \item{inf.mat}{Information matrix for the components of \code{delta}.}
   \item{var.coef}{Variance covariance matrix for \code{coef}.}
}

\author{
Robert E. Fay, Mamadou Diallo
}
