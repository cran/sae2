% File ~/sae2/man/eblupDyn.Rd
\name{eblupDyn}
\alias{eblupDyn}
\alias{eblupRY}
\title{EBLUP Fit of the Dynamic and Rao-Yu Time Series Models}
\description{
  Functions for producing EBLUP small area estimates of the dynamic or
  Rao-Yu time series models through either ML or REML estimation of the variance
  components. The functions can fit univariate or multivariate models.
}
\usage{
eblupDyn(formula,  D,  TI,  vardir,  method = c("REML", "ML"),
         MAXITER = 1000,  PRECISION = .1e-05,  data, 
         max.rho = NULL,  dampening = 1, \dots) 
         
         
eblupRY(formula,  D,  TI,  vardir,  method = c("REML", "ML"),
         MAXITER = 1000,  PRECISION = .1e-05,  data, 
         max.rho = .98, dampening = 0.9, \dots)
}
\arguments{
  \item{formula}{For a univariate model, a \code{formula} for the linear
    regression relationship between the dependent variable and the
    independent variable(s). The variables included in formula must have
    length equal to \code{D*TI} and be sorted in ascending order by time 
    within each domain.
    
    
    For a multivariate model, a list of formulas, one for each 
    dependent variable. The number of dependent variables, \code{NV}, is 
    determined from the length of the list. The dependent variables included 
    in the formulas must each have length equal to \code{D*TI} and be sorted
    in ascending order by time within each component within each domain,
    which is an extension of the sorting requirement for the univariate 
    model. Further details of the model specification are given under 
    Details.
    }
  \item{D}{The total number of domains.}
  \item{TI}{The number of time instances, typically years (constant for all 
    domains).}
  \item{vardir}{For the univariate model, the sampling covariance matrix 
    for the direct estimates of the \code{D*TI} elements of the dependent 
    variable. The covariance matrix should be in the form of a square
    matrix with \code{D*TI} rows and columns. Non-zero covariances between 
    domains are not allowed, so the matrix must have a block diagonal form 
    with \code{D} blocks, each of which is a square matrix with \code{TI} 
    rows and columns. Note that within domain, non-zero covariances are
    allowed over time. 
    
    Alternatively, \code{vardir} can be a list of 
    \code{D} covariance matrices, each with \code{TI} rows and columns.  

    
    For the multivariate model, the square covariance matrix for the 
    \code{D*NV*TI} elements of the dependent variables. The matrix must be
    in the form of a square matrix with \code{D*NV*TI} rows and columns. The 
    variances and covariances should be in the sort order of time within
    dependent variable within domain. Non-zero covariances between domains 
    are not allowed, but non-zero covariances may be present across time 
    and between components. 
    
    Alternatively, \code{vardir} can be a list of
    \code{D} covariance matrices, each with \code{NV*TI} rows and columns.}
  \item{method}{Whether restricted maximum likelihood \code{REML} or 
    maximum likelihood \code{ML} should be used.}
  \item{MAXITER}{The maximum number of iterations allowed for the 
    Fisher-scoring algorithm.}
  \item{PRECISION}{The convergence tolerance limit for the 
    Fisher-scoring algorithm.}
  \item{data}{An optional data frame containing the variables named in 
    \code{formula}. By default the variables are taken from the 
     environment from which \code{eblupDyn} is called. Because 
     \code{vardir} will be of a different size than the variables
     in \code{formula}, \code{data} will not be searched for
     \code{vardir}.}
  \item{max.rho}{If not \code{NULL}, the maximum value allowed for 
     \code{rho}. Note the different defaults for \code{eblupDyn} and
     \code{eblupRY}.}
  \item{dampening}{A multiplier of the computed update to parameters
     after iteration 5. A value less than 1 may slow the iterations but 
     lessens the chance of overshooting the optimum choice. The default 
     values were determined experimentally, but may be modified.}
  \item{\dots}{Other parameters passed to \code{dynRYfit} that affect
       convergence, provide starting values, or request specific results. 
       The exceptions are \code{y}, \code{X}, \code{NV}, \code{M}, 
       \code{ncolx}, and \code{model}, which will be set by either
       \code{eblupDyn} or \code{eblupRY}.}
}
\details{
A typical model has the form \code{response ~ terms} where \code{response}
is the (numeric) response vector and \code{terms} is a series of terms that
specifies a linear predictor for response. 


A formula has an implied intercept term. To remove this use either 
y ~ x - 1 or y ~ 0 + x. See \code{\link{formula}} for more details of 
allowed formulae.


\code{eblupDyn} and \code{eblupRY} parse \code{formula} by calling core
R functions to determine \code{X}, then calling \code{\link{dynRYfit}}.
As a last step, \code{eblupDyn} or \code{eblupRY} finalize the list that 
they return.


The additional parameters passed to \code{\link{dynRYfit}} may  
include \code{contrast.matrix}, which  specifies linear combinations 
of estimates within domains, such as the sum over dependent variables 
or moving averages across time. Corresponding MSE estimates are provided 
for the contrasts.


The argument \code{ids} accepts a data frame with \code{D}
rows of domain identifiers. These ids are returned in the list from
\code{eblupDyn} or \code{eblupRY}. 


If \code{iter.history} is set to \code{TRUE}, the returned object will 
include additional items with values of statistics at each step of the 
iteration; see \code{\link{dynRYfit}} for details on \code{delta.hist}, 
\code{llikelihood.hist}, \code{adj.hist}, \code{inf.mat.hist},
\code{s.hist}, \code{ix.hist}, \code{adj.factor.hist}, and 
\code{warning.hist}. The default action
is to include the history only if the iterations fail, in which case
the history might suggest what went wrong. In the case of convergence,
the history is usually not of interest, in which case omitting it
reduces the size of the returned object.


MSE estimation for REML for both the Rao-Yu and dynamic models follows
the results summarized in Rao and Molina (2015, pp. 98-111). The MSE estimates 
incorporate g1, g2, and g3 terms. Our simulations show that the REML estimates
have somewhat smaller MSEs than the ML estimates, but this is not reflected in 
the comparison of the estimated MSEs returned by the functions. The MSE 
estimates under REML perform quite well on average. The MSE estimates for ML 
use the same estimator as for REML, but they are modest underestimates of the 
true MSE in the same simulations.

}
\value{
   \item{eblup}{In the univariate case, a vector of length \code{D*TI} with 
   the eblup estimates. In the multivariate case, a data frame of D*TI rows 
   and NV columns.} 
   \item{fit}{A list summarizing the fit of the model with the following:
     \itemize{
       \item \code{model:} form of the model: T - Dynamic or Rao-Yu; REML
            or ML.
       \item \code{covergence:} a logical value indicating whether the 
                 convergence criterion was met.
       \item \code{iterations:} number of iterations performed by the 
                Fisher-scoring algorithm.
       \item \code{estcoef:} a data frame with the estimated model 
           coefficients (\code{beta}) in the first column , 
           their asymptotic standard errors (\code{std.error}) in the 
           second column, the t statistics (\code{tvalue}) in the third column, 
           and the p-values (\code{pvalue}) of the significance of each 
           coefficient in last column.
       \item \code{estvarcomp:}  a data frame with the estimated values 
           of the variances and correlation coefficients in the first column
           (\code{estimate}) and their asymptotic standard errors in the 
           second column (\code{std.error}).
       \item \code{goodness:} the log-likelihood and, if REML, the restricted 
           log-likelihood.
       }
     }
   \item{parm}{A labelled vector with the estimated variance components, 
         correlations, and number of iterations.}
   \item{coef}{A labelled vector of coefficients of the model or models.}
   \item{ids}{A data frame with \code{D} rows and one or more columns of 
         numeric or character domain identifiers.}
   \item{delta}{An ordered vector of the variance components, which may be 
         used as starting values for additional iterations, see 
         \code{dynRYfit}.}
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
   \item{model}{The formula or list of formulas implemented.}
}

\references{
- Fay, R.E. and Herriot, R.A. (1979). Estimation of income from small places: An application of James-Stein procedures to census data. Journal of the American Statistical Association 74, 269-277.

- Fay, R.E., Planty, M. and Diallo, M.S. (2013). Small area estimates from the National Crime Victimization Survey. Proceedings of the Joint Statistical Meetings. American Statistical Association, pp. 1544-1557.

- Rao, J.N.K. and Molina, I. (2015). Small Area Estimation, 2nd ed. Wiley, Hoboken, NJ.

- Rao, J.N.K. and Yu, M. (1994). Small area estimation by combining time series and cross-sectional data. Canadian Journal of Statistics 22, 511-528.

}
\author{
Robert E. Fay, Mamadou Diallo
}
\examples{
D <- 20 # number of domains
TI <- 5 # number of years
set.seed(1)
data <- data.frame(Y= mvrnormSeries(D=D, TI=TI, rho.dyn=.9, sigma.v.dyn=1, 
   sigma.u.dyn=.19, sigma.e=diag(5)), X=rep(1:TI, times=D))
result.dyn  <- eblupDyn(Y ~ X, D, TI, vardir = diag(100), data=data)
result.dyn$fit

require(sae)
data(spacetime)      # Load data set from sae package
data(spacetimeprox)  # Load proximity matrix 

D <- nrow(spacetimeprox)            # number of domains
TI <- length(unique(spacetime$Time)) # number of time instants
# Fit model ST with AR(1) time effects for each domain
resultST <- eblupSTFH(Y ~ X1 + X2, D, TI, Var, spacetimeprox,
                      data=spacetime)
resultT  <- eblupDyn(Y ~ X1 + X2, D, TI, vardir = diag(spacetime$Var),
                      data=spacetime, ids=spacetime$Area)
resultT.RY  <- eblupRY(Y ~ X1 + X2, D, TI, vardir = diag(spacetime$Var),
                      data=spacetime, ids=spacetime$Area)
resultST$fit
resultT$fit
resultT.RY$fit
rowsT <- seq(TI, TI*D, by=TI)
data.frame(Domain=spacetime$Area[rowsT], Y=spacetime$Y[rowsT], 
              EBLUP_ST=resultST$eblup[rowsT],
              EBLUB_Dyn=resultT$eblup[rowsT],
              EBLUP_RY=resultT.RY$eblup[rowsT])
}

