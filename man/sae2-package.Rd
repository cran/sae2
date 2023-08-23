\name{sae2-package}
\alias{sae2-package}
\alias{sae2}
\docType{package}
\title{
Small Area Estimation: Time-Series Models.
}
\description{
Time series are-level models for small area estimation. 
      The package supplements the functionality of the package \code{sae}. 
      Specifically, it includes EBLUP fitting of the original 
      Rao-Yu model, which did not have a spatial component. It also offers 
      a modified ("dynamic") version of the Rao-Yu model, replacing the 
      assumption of stationarity. Both univariate and multivariate applications 
      are supported. Of particular note is the allowance for covariance of
      the area-level sample estimates over time, as encountered in
      rotating panel designs such as the U.S. National Crime Victimization
      Survey or present in a time-series of 5-year estimates from the 
      American Community Survey.
}
\details{
\tabular{ll}{
Package: \tab sae2\cr
Type: \tab Package\cr
Version: \tab 1.2-1\cr
Date: \tab 2023-08-22\cr
License: \tab GPL-2\cr
}
The package provides two primary functions, \code{eblupRY} and \code{eblupDyn}, 
to fit non-spatial time-series small area models to area-level data. The 
function \code{mvrnormSeries} provides simulated data under either model. 
Functions \code{geo_ratios} and \code{vcovgen} can assist in preparing the
input.
}
\author{
Robert E. Fay, Mamadou S. Diallo

Maintainer: Robert E. Fay <bobfay@hotmail.com>

}
\references{
- Fay, R.E. and Herriot, R.A. (1979). Estimation of income from small places: An application of James-Stein procedures to census data. Journal of the American Statistical Association 74, 269-277.

- Fay, R.E., Planty, M. and Diallo, M.S. (2013). Small area estimates from the National Crime Victimization Survey. Proceedings of the Joint Statistical Meetings. American Statistical Association, pp. 1544-1557.

- Rao, J.N.K. and Molina, I. (2015). Small Area Estimation, 2nd ed. Wiley, Hoboken, NJ.

- Rao, J.N.K. and Yu, M. (1994). Small area estimation by combining time series and cross-sectional data. Canadian Journal of Statistics 22, 511-528.

}
\seealso{
\code{\link[sae:sae-package]{sae}}
}
