# - gamma_u and its derivative for the Rao-Yu model 

Gamma_u_ry <- function(Rho, TI) {
   if (length(TI) > 1) {
      nt <- tail(TI, n=1)
   } else {
      nt <- TI
   }
   Gamma_b <- matrix(rep(1:nt,times=nt),nrow=nt, ncol=nt)
   Gamma_s <- abs(Gamma_b - t(Gamma_b))
   if(Rho > 0) {
      Gamma_u <- Rho^(Gamma_s)/(1-Rho^2)
      d.Gamma_u <- ((Gamma_s * Rho^(Gamma_s) * (1-Rho^2)/Rho) +
                   2 * Rho * Rho^(Gamma_s))/((1-Rho^2)^2)
   } else {
      Gamma_u <- diag(nt)
      d.Gamma_u <- matrix(0, nrow=nt, ncol=nt)
      for (i in 2:nt) {
         d.Gamma_u[i-1, i] < -1
         d.Gamma_u[i, i-1] < -1 
      }
   }
   if (length(TI) > 1) {
      Gamma_u <- Gamma_u[TI, TI]
      d.Gamma_u <- d.Gamma_u[TI, TI]
   }
   return(list(Gamma_u= Gamma_u, d.Gamma_u = d.Gamma_u)) 
} 
