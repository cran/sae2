% File ~/sae2/inst/NEWS.Rd
\name{sae2NEWS}
\title{sae2 NEWS}
\section{Changes in version 1.2-2 (2025-08-26)}{
     \itemize{
        \item  The argument TI in the functions eblupDyn, eblupRY, and dynRYfit has 
        a new downward-compatible generalization. TI of length 1 functions as 
        previously. TI with a length greater than 1 indexes the years in the series to be 
        included in the model, supporting a model with some breaks in the series.
        
        \item The parameter y.scale, if specified and not equal to 1, previously caused 
        incorrect estimation of the mse.
        
        \item Minor changes have been made to the default convergence criterion,
        and more detailed error messages are provided from 1.2-1. 
     }
}