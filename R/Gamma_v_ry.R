# - gamma_v and its derivative for the Rao-Yu model 

Gamma_v_ry <- function(Rho, TI) {
   if (length(TI) > 1) {
      nt <- tail(TI, n=1)
   } else {
      nt <- TI
   }
   Gamma_v <- matrix(1, nrow=nt, ncol=nt) 
   d.Gamma_v <- matrix(0, nrow=nt, ncol=nt)
   if (length(TI) > 1) {
      Gamma_v <- Gamma_v[TI, TI]
      d.Gamma_v <- d.Gamma_v[TI, TI]
   }
   return(list(Gamma_v= Gamma_v, d.Gamma_v = d.Gamma_v)) 
} 
