\name{sae2-package}
\alias{sae2-package}
\alias{sae2}
\docType{package}
\title{
Small Area Estimation: Time-series Models.
}
\description{
Time series models to supplement the functionality of the package \code{sae}, 
specifically EBLUP fitting of the original Rao-Yu model and also of a modified
("dynamic") version. Both univariate and multivariate applications are 
supported. The models are non-spatial, in contrast to the spatial models in 
\code{sae}.}
\details{
\tabular{ll}{
Package: \tab sae2\cr
Type: \tab Package\cr
Version: \tab 0.1-1\cr
Date: \tab 2015-01-26\cr
License: \tab GPL-2\cr
}
The package provides two primary functions, \code{eblupRY} and \code{eblupDyn}, 
to fit non-spatial time-series small area models to area-level data. The 
function \code{mvrnormSeries} provides simulated data under either model.
}
\author{
Robert E. Fay, Mamadou S. Diallo

Maintainer: Robert E. Fay <bobfay@hotmail.com>

}
\references{
- Fay, R.E. and Herriot, R.A. (1979). Estimation of income from small places: An application of James-Stein procedures to census data. Journal of the American Statistical Association 74, 269-277.

- Fay, R.E., Planty, M. and Diallo, M.S. (2013). Small area estimates from the National Crime Victimization Survey. Proceedings of the Joint Statistical Meetings. American Statistical Association, pp. 1544-1557.

- Rao, J.N.K. (2003). Small Area Estimation. Wiley, New York.

- Rao, J.N.K. and Yu, M. (1994). Small area estimation by combining time series and cross-sectional data. Canadian Journal of Statistics 22, 511-528.

}
\seealso{
\code{\link[sae:sae-package]{sae}}
}
