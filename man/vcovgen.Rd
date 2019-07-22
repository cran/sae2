% File c:/sae2/man/vcovgen.Rd
\name{vcovgen}
\alias{vcovgen}
\title{Estimate variance/covariance matrices using linear substitutes 
}
\description{The function computes estimates of variance/covariance 
       matrices for sampling with replacement at the first stage
       from a stratified sample, using the linear substitutes
       produced by \code{\link{geo_ratios}}. The function produces a 
       list of estimated covariance matrices corresponding the 
       geographic areas identified by \code{geocode}. 
}
\usage{
vcovgen(linear.subs, year.list, geocode, designvars)

}
\arguments{
  \item{linear.subs}{A data frame with linear substitutes for each
                     Y variable, \code{Year}, and the variables with
                     the names given by \code{geocode} and 
                     \code{designvars}.
}
  \item{year.list}{A character or numeric vector of the years to 
                   include.
}
  \item{geocode}{A character variable with the name of the geographic 
                 variable for which separate estimates are of interest.
}
  \item{designvars}{A character vector of length 2, naming the strata 
                and first-stage units.                    
}
}
\details{
The function reformats the second element of the list output from
\code{geo_ratios} in order to estimate covariances across time for 
each separate geographic area identified by the variable with the
name indicated by \code{geocode}. It then calls functions from the
\pkg{survey} package to estimate covariances.

See the example in \code{\link{geo_ratios}}.
}
\value{A list of estimated covariance matrices, sorted by the geographic
names. The list may be used by \code{eblupDyn} or \code{eblupRY}, although
many applications require the estimated variance/covariance matrices to
be smoothed.
}
\author{Robert E. Fay
}
