# - gamma_u and its derivative for the dynamic model 

Gamma_u_f <- function(Rho, TI) {
   if (length(TI) > 1) {
      nt <- tail(TI, n=1)
   } else {
      nt <- TI
   }
   Gamma_b <- matrix(rep(1:nt, times=nt), nrow=nt, ncol=nt) 
   Gamma <- Rho^(abs(Gamma_b - t(Gamma_b)))
   diag_u <- rep(0, times=nt)
   d.diag_u <- diag_u
   for (i in 2:nt) {
       diag_u[i] <- diag_u[i-1] + Rho^(2*i-4)
       d.diag_u[i] <- d.diag_u[i-1] + (2*i-4)*Rho^(2*i-4)/Rho
       }
   Gamma_u <- diag(diag_u)
   d.Gamma_u <- diag(d.diag_u)
   d.Gamma <- abs(Gamma_b - t(Gamma_b))*Gamma/Rho
   for (i in 2:(nt-1)) {
       for (j in (i+1):nt) {
          Gamma_u[i,j] <- diag_u[i] * Gamma[i,j]
          Gamma_u[j,i] <- Gamma_u[i,j]
          d.Gamma_u[i,j] <- d.diag_u[i] * Gamma[i,j] + 
                            diag_u[i] * d.Gamma[i,j]
          d.Gamma_u[j,i] <- d.Gamma_u[i,j]
        }
     }
   if (length(TI) > 1) {
      Gamma_u <- Gamma_u[TI, TI]
      d.Gamma_u <- d.Gamma_u[TI, TI]
   }
   return(list(Gamma_u= Gamma_u, d.Gamma_u = d.Gamma_u)) 
}

