# - gamma_u and its derivative for the dynamic model 

Gamma_u_f <- function(Rho,TI) {
  Gamma_b <- matrix(rep(1:TI,times=TI),nrow=TI, ncol=TI) 
  Gamma <- Rho^(abs(Gamma_b-t(Gamma_b)))
  diag_u<-rep(0,times=TI)
  d.diag_u<-diag_u
  for (i in 2:TI) {
      diag_u[i] <-diag_u[i-1]+Rho^(2*i-4)
      d.diag_u[i]<-d.diag_u[i-1]+(2*i-4)*Rho^(2*i-4)/Rho
      }
  Gamma_u <- diag(diag_u)
  d.Gamma_u <- diag(d.diag_u)
  d.Gamma <- abs(Gamma_b-t(Gamma_b))*Gamma/Rho
  for (i in 2:(TI-1)) {
      for (j in (i+1):TI) {
         Gamma_u[i,j] <- diag_u[i] * Gamma[i,j]
         Gamma_u[j,i] <- Gamma_u[i,j]
         d.Gamma_u[i,j] <- d.diag_u[i] * Gamma[i,j] + 
                           diag_u[i] * d.Gamma[i,j]
         d.Gamma_u[j,i] <- d.Gamma_u[i,j]
         }
      }
  list(Gamma_u= Gamma_u, d.Gamma_u = d.Gamma_u) 
  }

