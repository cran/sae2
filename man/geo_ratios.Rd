% File c:/sae2/man/geo_ratios.Rd
\name{geo_ratios}
\alias{geo_ratios}
\title{Compute rates or ratios for a set of geographic entities over a set 
       of years
}
\description{The function computes rates or ratios by a geographic code and 
       the variable \code{Year}. If \code{designvars} is specified, the function 
       also returns a data frame with linear substitutes to compute Taylor 
       series variances.
}
\usage{
geo_ratios(data,  geocode,  numerators,  denominators,  geonames, 
           new.names,  designvars)

}
\arguments{
  \item{data}{A data frame with the required variables, including a variable 
           named \code{Year}.
}
  \item{geocode}{A character variable with the name of the geographic 
                 variable for which separate estimates are of interest.
}
  \item{numerators}{A character vector listing the names in \code{data} of the
                  numerators of the ratios.
}
  \item{denominators}{A character vector listing the names in \code{data} of 
                 the denominators of the ratios. If a single value is given,
                 it  will be used for all of the ratios.
}
  \item{geonames}{An optional data frame containing \code{geocode} and one or
                  more geographic variables, such as names, that will be merged
                  into the results. There should be only one row for each
                  value of \code{geocode}.
}
      
  \item{new.names}{An optional character vector of the same length as 
                  \code{numerators} naming the resulting ratios. 
                  If \code{new.names} is not specified, the output ratios
                  will have the same names as \code{numerators}.
}
  \item{designvars}{Optional. If given, a character vector naming one or more 
                  survey design variables in \code{data} to use in forming 
                  linear substitutes for variance calculation of the ratios 
                  cross-classified by \code{Year} and the variable named by 
                  \code{geocode}. The vector should not include the 
                  \code{geocode} variable or \code{Year}. 
}
}
\details{
For programming simplicity, the function enforces the requirement that names 
should not be repeated in either \code{numerators} or \code{new.names}. Names 
may be repeated in \code{denominators}.


Rather than a typical survey file, the function expects the data frame 
\code{data} to contain weighted estimates for each analytic variable. As a 
simple example, to find the variance of the weighted mean of \code{y} with
weights \code{w}, \code{data} should contain \code{w} and 
\code{y * w}. For convenience, the weighted estimates can still be 
assigned their original names in \code{data}, such as \code{y}. In this 
case,

   \code{numerators = y, denominators=w}

would create the appropriate linear substitutes for the variance of the 
weighted mean.

   
This design of the function allows complex possibilities, such as estimating
the variance of a rate where the numerator is based on one weight and the 
denominator is based on another. For example, estimation for the National 
Crime Victimization Survey requires this capability.
}
\value{If \code{designvars} is not specified, a named list with one element,
                  a data frame containing the ratios sorted by \code{geocode}
                  and \code{Year}. 
                  
                  If \code{designvars} is specified, a second element is added 
                  to the list, a data frame giving the totals of the linear 
                  substitutes by \code{Year}, \code{geocode}, and 
                  \code{designvars}. The elements of the list are named 
                  \code{estimates} and \code{linear.subs}.               
}

\references{
- Woodruff, R.S. (1971). A simple method for approximating the variance of a complex estimate. Journal of the American Statistical Association 66, 411-414.

}

\author{Robert E. Fay
}
\examples{
require(survey)
require(MASS)
D <- 20 # number of domains
T <- 5 # number of years
samp <- 16 # number of sample cases per domain
set.seed(1)
# use conditional.mean=TRUE to generate true small area values
# without sampling error
Y.list <- mvrnormSeries(D=D, T=T, rho.dyn=.9, sigma.v.dyn=1, 
   sigma.u.dyn=.19, sigma.e=diag(5), conditional.mean=TRUE)
# generate sampling errors
e <- rnorm(samp * T * D, mean=0, sd=4)
Y <- Y.list[[2]] + tapply(e, rep(1:100, each=16), mean) 
data <- data.frame(Y=Y, X=rep(1:T, times=D))
# model fit with the true sampling variances
result.dyn  <- eblupDyn(Y ~ X, D, T, vardir = diag(100), data=data)
# individual level observations consistent with Y
Y2 <- rep(Y.list[[2]], each=16) + e 
data2 <- data.frame(Y=Y2, X=rep(rep(1:T, each=samp), times=D), 
                   Year=rep(rep(1:T, each=samp), times=D), 
                   weight=rep(1, times=samp*T*D),
                   d=rep(1:D, each=samp*T), 
                   strata=rep(1:(D*T), each=samp),
                   ids=1:(D*T*samp))
# geo_ratios with designvars specified
geo.results <- geo_ratios(data2, geocode="d", numerators="Y",
                          denominators="weight",
                          designvars=c("strata", "ids"))
# illustrative check                          
max(abs(geo.results[[1]]$Y - Y))
vcov.list <- vcovgen(geo.results[[2]], year.list=1:5, geocode="d", 
      designvars=c("strata", "ids"))
vcov.list[[1]]
# model fitted with directly estimated variance-covariances
result2.dyn <- eblupDyn(Y ~ X, D, T, vardir=vcov.list, data=data)
cor(result.dyn$eblup, result2.dyn$eblup)
}
\seealso{\code{vcovgen}}