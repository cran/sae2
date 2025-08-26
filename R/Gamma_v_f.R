# - gamma_v and its derivative for the dynamic model 

Gamma_v_f <- function(Rho, TI) {
   if (length(TI) > 1) {
      nt <- tail(TI, n=1)
   } else {
      nt <- TI
   }
   rho_power <- Rho^(0:(nt-1))
   rho_d <- c(0,(Rho^(0:(nt-2)))*(1:(nt-1)) )
   Gamma_v= rho_power %*% t(rho_power)  
   d.Gamma_v = rho_power %*% t(rho_d) +
                  rho_d %*% t(rho_power)
   if (length(TI) > 1) {
      Gamma_v <- Gamma_v[TI, TI]
      d.Gamma_v <- d.Gamma_v[TI, TI]
   }
   return(list(Gamma_v= Gamma_v,
        d.Gamma_v = d.Gamma_v))
 }
