% File ~/sae2/man/remldyn.Rd
\name{dynRYfit}
\alias{dynRYfit}
\title{Internal fitting function for Dynamic and Rao-Yu models}
\description{
  Function designed to be called by either \code{eblupDyn} or \code{eblupRY}
  to produce EBLUP small area estimates of the dynamic or Rao-Yu time 
  series models through either ML or REML estimation of the variance 
  components. For completeness, the function is documented here, but
  users are encouraged to base applications on calls to the more
  convenient \code{eblupDyn} or \code{eblupRY}.
}
\usage{
dynRYfit(y,  X,  M,  T,  NV=1,  vcov_e,  maxiter=100, 
         iter.tol=.1e-5,  ncolx=NULL,  sig2_u=1,  sig2_v=1, 
         rho=.8,  rho_u =.4,  delta=NULL,  rho.fixed=NULL, 
         y.include=NULL, ids=NULL,  contrast.matrix=NULL,  
         baby.steps=TRUE,  dampening=NULL,  iter.history=NULL,  
         sig2.min.factor=.0001,  max.rho_u=.98,  max.rho=NULL, 
         tol=.Machine$double.eps, y.rescale=NULL, 
         llike.only=FALSE,  method=c("REML", "ML"), 
         model=c("dyn", "RY"))
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
    \code{M*NV*T} elements of the dependent variables. The matrix should be 
    in the form of a square matrix with \code{M*NV*T} rows and columns. Time 
    should vary within variable, which should vary within domain.
    Non-zero covariances between domains are not allowed, but non-zero 
    covariances may be present across time and between variables.}
  \item{maxiter}{The maximum number of iterations allowed for the 
    Fisher-scoring algorithm, with a default value of 100.}
  \item{iter.tol}{The convergence tolerance limit for the 
    Fisher-scoring algorithm, with a default value of .000001.}
  \item{ncolx}{For a univariate model, the number of columns of X. For a 
    multivariate model, a vector of length NV must be specified giving the
    number of columns of X used for each dependent variable.}
  \item{sig2_u}{An initial starting value or values for the variance of the
    random increments.}
  \item{sig2_v}{An initial starting value or values for a domain level
    random effect. In the Rao-Yu model, the random effect is constant over
    time, whereas in the dynamic model it is an initial effect subject to 
    dampening over time.}
  \item{rho}{The correlation across time. This correlation is assumed to
    be the same for the dependent variables in the multivariate model.}
  \item{rho_u}{For \code{NV > 1} only, the (NV*(NV-1))/2 starting values for
    the correlations between the random effects of the different dependent
    variables. If a single value is given, it will be used for the 
    (NV*(NV-1))/2 components. The sort order corresponds to a lower triangle
    of the covariance matrix.}
  \item{delta}{The random effect components in the preferred internal order.
    Specification of \code{delta} will override any specification of 
    \code{sig2_u}, \code{sig2_v}, \code{rho}, or \code{rho_u}. In the 
    univariate case, delta should contain \code{sig2_u}, \code{sig2_v}, and 
    \code{rho}. In the bivariate case, \code{delta} should be of length 6 with
    \code{sig2_u} and \code{sig2_v}, each of length 2, \code{rho}, and 
    \code{rho_u}. Similarly, for 3 dependent variables, the length of 
    \code{delta} is 10, with the last 3 elements \code{rho_u} in lower
    triangular form.}
  \item{rho.fixed}{If TRUE, the value of \code{rho} imbedded in \code{delta},
    if specified, or else given by \code{rho} will remain fixed during the
    iterations. Among other features, this allows the likelihood function
    for trial values of \code{rho} to be computed at the maximum over the
    other random effect parameters.}
  \item{y.include}{If specified, vector of length \code{M} to indicate 
    which domains to include in the estimation, with 1 signalling inclusion
    and 0 exclusion. Estimates for the excluded domains will be based on
    the fixed effects model only.}
  \item{ids}{A data frame with \code{M} rows giving ids for each of the
             domains. The data frame is copied into the returned object.}
  \item{contrast.matrix}{A matrix of coefficients of contrasts. The matrix
             must have \code{T*NV} rows, but it can contain an arbitrary number
             of columns. Within each domain, the coefficients are applied
             to the \code{T*NV} EBLUP estimates.}
  \item{baby.steps}{Unless specified as FALSE, the first five iterations of the
             Fisher scoring algorithm are dampened by factors of 
             \code{c(.0625, .125, .25, .5, .75)}. These heuristically derived 
             factors appear to lessen drastic overshooting of the true maximum
             in the initial iterations.}
  \item{dampening}{A factor used to dampen the changes in the random
             effect parameters. Unlike \code{baby.steps}, its effect persists
             during all of the iterations until convergence. Note that the 
             "factory setting" of this parameter is 1 for the dynamic model but
             .9 for the Rao-Yu model.}
  \item{iter.history}{If TRUE, key values are saved during each  
             iteration and included as additional items, described below,
             in the returned list: \code{delta.hist}, \code{llikelihood.hist},
             \code{adj.hist}, \code{s.hist}, \code{ix.hist}, 
             \code{adj.factor.hist}, and \code{warning.hist}.
             
             If \code{iter.history} is not specified, these items will be 
             returned only if the calculations are begun but not successfully 
             completed.}
  \item{sig2.min.factor}{A factor to multiply the minimum direct variance to
         use as a minimum value for any of the variance components. The 
         iterations will be constrained not to go below the resulting bounds.}
  \item{max.rho_u}{A maximum allowed value for the estimated \code{rho_u}. The
         default value is .98.}
  \item{max.rho}{A maximum allowed value for \code{rho}. By default, 
         \code{\link{eblupRY}} sets this value to .98.}
  \item{tol}{A tolerance value used by matrix routines to prevent numerical
         instability. The value may be set to a lower value to encourage
         covergence, but appropriate caution should be applied.}
  \item{y.rescale}{In the univariate case, a scaler multiplier for all of
         the \code{y} values. If the \code{y} values are either too small or 
         too large, the information matrix may become unstable. Setting this
         value to 1 has no effect; setting it to 10 or 100 rescales very 
         small \code{y} values to a more appropriate range. Similarly, 
         positive values less than 1 may be used to rescale large \code{y} 
         values. The effect of rescaling is removed before normal return 
         from the function, within the limits of normal precision.
         
         In the multivariate case, \code{y.rescale} may be a vector of 
         length \code{NV} to rescale each component separately.}
  \item{llike.only}{Compute the log-likelihood (ML) or restricted 
         log-likelihood (REML) without further iteration, typically from
         values specified by \code{delta}.}
  \item{method}{Use restricted maximum likelihood (\code{"REML"}) or 
         maximum likelihood (\code{"ML"}).}
  \item{model}{Dynamic (\code{"dyn"}) or Rao-Yu (\code{"RY"}).}
}
\details{
   Many of arguments can be used to control the iterations if the 
   defaults lead to convergence difficulties. 
   
   \code{llike.only} in combination with \code{delta} permits a 
    point-by-point investigation of the likelihood surface. 
   
   The primary functions \code{eblupDyn} and \code{eblupRY} determine \code{X}, 
   \code{NV}, \code{colx}, and \code{model}, but the remaining parameters can 
   be passed to \code{dynRYfit} through \code{eblupDyn} or \code{eblupRY}.
}
\value{
   \item{eblup}{In the univariate case, a vector of length \code{M*T} with 
   the eblup estimates in the same sort order as \code{y}. In the multivariate
   case, a matrix of M*T rows and NV columns.}
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
   \item{coef}{A labelled vector of coefficients for the fixed effects of 
         the model or models.}
   \item{ids}{A data frame with \code{D} rows and one or more columns of 
         numeric or character domain identifiers.}
   \item{delta}{An ordered vector of the variance components (see above). 
         It may be used as starting values for additional iterations.}
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
   \item{delta.hist}{Values of \code{delta} at each iteration.}
   \item{llikelihood.hist}{Values of the log-likelihood (ML) or
       restricted log-likelihood (REML) at each iteration.}
   \item{adj.hist}{Number of cycles in the internal loop to determine 
       \code{adj.factor} within each iteration.}
   \item{inf.mat.hist}{Values of \code{inf.mat} at each iteration.}
   \item{s.hist}{Vector to be multiplied by the inverse information matrix
       to determine the change in the parameters.}
   \item{ix.hist}{List of parameters eligible for change at each iteration.
       Parameters with estimated changes out of bounds will not be 
       eligible.}
   \item{adj.factor.hist}{Adjustment to the vector change in the parameters 
       at each iteration.}
   \item{warning.hist}{A 4-row matrix of warnings at each iteration, where
       warning 1 is set to 1 for the iteration if the algorithm has not found 
       an increase in the restricted log likelihood or log likelihood, 
       warning 2 is set to 1 if the maximum number of iterations is reached, 
       warning 3 is set to 1 if the estimated variance-covariance matrix 
       becomes singular, and warning 4 is set to 1 if the coefficients of the
       fixed effects cannot be estimated.}
}

\author{
Robert E. Fay, Mamadou Diallo
}
