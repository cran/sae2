# dynRYfit - maximum likelihood or restricted maximum likelihood estimation 
# of the multivariate dynamic or Rao-Yu model

dynRYfit <- function(y, X, M, TI, NV=1, vcov_e, maxiter=100, iter.tol=.1e-5,
         ncolx=NULL, sig2_u=1, sig2_v=1, rho=.8, rho_u =.4, delta=NULL,  
         rho.fixed=NULL, y.include=NULL, ids=NULL, contrast.matrix=NULL,  
         baby.steps=TRUE, dampening=NULL, iter.history=NULL,
         sig2.min.factor=.0001, max.rho_u=.98, max.rho=NULL, 
         tol=.Machine$double.eps, y.rescale=NULL, llike.only=FALSE, 
         method=c("REML", "ML"), model=c("dyn", "RY")) {
         
   if (length(model) > 1) model <- model[1]
   if (!model %in% c("dyn", "RY"))
      stop(" model=\"", model, "\" must be \"dyn\", or \"RY\"")
   if (length(method) > 1) method <- method[1]
   if (!method %in% c("REML", "ML", "MLE")) 
      stop(" method=\"", method, "\" must be \"REML\", or \"ML\"")
   if (model == "dyn") {
      if (method == "REML") {
         fit <- list(model="T: Dynamic, REML", convergence=FALSE)
      } else {
         fit <- list(model="T: Dynamic, ML", convergence=FALSE)
      }
   } else {
      if (method == "REML") {
         fit <- list(model="T: Rao-Yu, REML", convergence=FALSE)
      } else {
         fit <- list(model="T: Rao-Yu, ML", convergence=FALSE)
      }
   }   
   if(!is.null(y.rescale)) {
      if (any(y.rescale <= 0)) 
         stop("y.rescale must be positive")
      if (length(y.rescale) == 1) 
         y.rescale <- rep(y.rescale, NV)
      if (length(y.rescale) != NV) 
         stop("y.rescale must be of length 1 or NV")
      y.adjust <- rep(1, times=M*TI*NV)
      for (m in 1:M) {
         for (nv in 1:NV) {
            y.adjust[((m-1)*TI*NV + (nv-1)*TI + 1):((m-1)*TI*NV + nv*TI)] <-
                  y.rescale[nv]               
         }
      }    
      y <- y * y.adjust
      vcov_e <- y.adjust * t(vcov_e * y.adjust)
      use.y.rescale <- TRUE
   } else {
      use.y.rescale <- FALSE
      y.rescale <- rep(1, NV)
   }
   delta.adjust <- rep(1, 2*NV + 1 + (NV*(NV-1))/2) 
   delta.adjust[1:(2*NV)] <- c(y.rescale, y.rescale)
   istart <- 1
   b.adjust <- rep(0, dim(X)[2])
   for (nv in 1:NV) {
      b.adjust[istart:(istart+ncolx[nv]-1)] <- y.rescale[nv]
      istart <- istart + ncolx[nv]
   }
   if(is.null(delta)) {
      starting.delta <- FALSE 
   } else { 
      if (length(delta)!= 2*NV + 1 + (NV*(NV-1))/2) {
         starting.delta <- FALSE
         warning("delta of wrong length ignored", immediate.=TRUE) 
      } else { 
         starting.delta <- TRUE
         if(use.y.rescale)
            delta[1:(2*NV)] <- (delta.adjust[1:(2:NV)]^2) * delta[1:(2*NV)]
         sig2_u <- delta[1:NV]
         sig2_v <- delta[(NV+1):(2*NV)]
         if(is.null(rho.fixed)) {
            rho <- delta[2*NV+1]
         } else {  
            warning("Value of rho.fixed imposed on specified delta", 
               immediate.=TRUE)
            rho <- rho.fixed
         }    
         if(NV > 1) rho_u <- delta[(2*NV+2):(2*NV+1+(NV*(NV-1))/2)]
      }
   }

   if(!starting.delta){
      if(length(sig2_u)==1) sig2_u <- rep(sig2_u, NV)
      if(length(sig2_v)==1) sig2_v <- rep(sig2_v, NV)
      if(use.y.rescale) {
         sig2_u <- y.rescale^2 * sig2_u
         sig2_v <- y.rescale^2 * sig2_v
      }
      if(length(rho_u)==1) rho_u <- rep(rho_u, ((NV*(NV-1))/2))
      if(!is.null(rho.fixed)) rho <- rho.fixed
   }

   if(!is.null(y.include)) {
      if(length(y.include)!= M) stop("The length of y.include must match M")
      M.save <- M
      y.in <- which(y.include == 1)
      M <- length(y.in)
      if(M < M.save) {
         X.save <- X
         y.save <- y 
         vcov_e.save <- vcov_e
         index.obs <- rep(0, NV*M*TI)
         for (i in 1:M) {
            index.obs[(NV*TI*(i-1)+1):(NV*TI*i)] <-
                           c((NV*TI*(y.in[i]-1)+1):(NV*TI*y.in[i]))
         }
         if(use.y.rescale) {
            y.adjust.save <- y.adjust
            y.adjust <- y.adjust.save[index.obs]
         }
         y <- y.save[index.obs]
         X <- X.save[index.obs, ]
         vcov_e <- vcov_e.save[index.obs, index.obs]
      }
   } else {
      y.in <- 1:M
   }
   nonposeigen <- vector(mode = "integer", length=0)
   largecondition <- vector(mode = "integer", length=0)
   for (m in 1:M) {
      vcovx <- vcov_e[((m-1)*NV*TI+1):(m*NV*TI),
                      ((m-1)*NV*TI+1):(m*NV*TI)]
      values <- eigen(vcovx, symmetric=TRUE, only.values=TRUE)$values
      if (values[NV*TI] <= 0) {
         nonposeigen <- append(nonposeigen, y.in[m])
      } else if (values[1]/values[NV*TI] > .1e14) {
         largecondition <- append(largecondition, y.in[m]) 
      }
   }
   if (length(nonposeigen) > 0) {
      warning("Sampling variance has zero or negative eigenvalues", 
                    immediate.=TRUE)
      print(c("Areas", paste(nonposeigen)))
   }
   if (length(largecondition) > 0) {
      warning("Sampling variance has large condition number", immediate.=TRUE)
      print(c("Areas", paste(largecondition)))
   }
   if(is.null(contrast.matrix)) {
      mi.mat <- diag(NV*TI)
      n.contr <- NV*TI
      contrast.mse <- contrast.g1 <- contrast.g2 <- contrast.g3 <- 
          contrast.fixed.var <- NULL
   } else{
      if(is.vector(contrast.matrix)) contrast.matrix <- matrix(contrast.matrix,
                                          nrow=length(contrast.matrix), ncol=1)
      if(!is.matrix(contrast.matrix)) 
         stop("contrast.matrix must be a matrix or a vector")
      if(dim(contrast.matrix)[1]!=NV*TI) 
         stop("The number of rows for contrast.matrix must agree with NV*TI")
      n.contr <- NV*TI + dim(contrast.matrix)[2]
      mi.mat <- matrix(0, nrow=NV*TI, ncol=n.contr)
      mi.mat[, 1:(NV*TI)] <- diag(NV*TI)
      mi.mat[,(NV*TI+1):n.contr] <- contrast.matrix
      contrast.mse <- contrast.g1 <- contrast.g2 <- contrast.g3  <- 
         contrast.fixed.var <- matrix(0, nrow=M, ncol=dim(contrast.matrix)[2])
   }
   contrast.est <- contrast.fixed.est <- contrast.wt1 <- contrast.wt2 <- 
          matrix(0, nrow=M, ncol=n.contr)
     
   N <- NV*M*TI
   len.rho_u <- (NV*(NV-1))/2
   if(!starting.delta) {
      if (NV > 1){ 
         delta <- c(sig2_u=sig2_u[1:NV], sig2_v=sig2_v[1:NV], rho=rho, 
                    rho_u=rho_u[1:len.rho_u]) 
      } else {
         delta <- c(sig2_u=sig2_u[1:NV], sig2_v=sig2_v[1:NV], rho=rho)
      }
   } else {
      if (NV > 1) {
         names(delta) <- c(paste0("sig_u", 1:NV), paste0("sig_v", 1:NV),
                           "rho", paste0("rho_u", 1:len.rho_u))
      } else {
         names(delta) <- c(paste0("sig_u", 1:NV), paste0("sig_v", 1:NV),
                           "rho")
      }
   }
   len.delta <- 2*NV+1 + len.rho_u
   index.rho <- 2*NV+1
   ix.start <- 1:len.delta
   ix.len <- len.delta
   if(!is.null(rho.fixed)) {
      ix.start <- ix.start[-index.rho]
      ix.len <- ix.len - 1
   }
   min.vcov <- min(diag(vcov_e))
   delta.min <- rep(sig2.min.factor * min.vcov, 2*NV)
   delta.min[index.rho] <- 0
   if (NV > 1) delta.min[(2*NV+2):len.delta] <- 0
   delta.hist <- matrix(0, nrow=len.delta, ncol=maxiter)
   llikelihood.hist <- rep(0, times=maxiter)
   adj.hist <- rep(0, times=maxiter)
   adj.hist[1] <- 1
   inf.mat.hist <- array(0, dim=c(len.delta, len.delta, maxiter))
   s.hist <- matrix(0, nrow=len.delta, ncol=maxiter)
   adj.factor.hist <- rep(0, times=maxiter)
   ix.hist <- matrix(0, nrow=len.delta, ncol=maxiter)
   warning.hist <- matrix(0, nrow=4, ncol=maxiter)
   rho <- delta[index.rho]
   rho_u <- delta[(index.rho+1):len.delta]
   if (model == "dyn") {
      Gamma.list <- Gamma_u_f(rho, TI)
      Gamma_u <- Gamma.list[[1]]
      d.Gamma_u <- Gamma.list[[2]]
      Gamma.list <- Gamma_v_f(rho, TI)
      Gamma_v <- Gamma.list[[1]]
      d.Gamma_v <- Gamma.list[[2]]
   } else {
      Gamma_b <- matrix(rep(1:TI,times=TI),nrow=TI, ncol=TI)
      Gamma_s <- abs(Gamma_b-t(Gamma_b))
      if(rho > 0) {
         Gamma_u <- rho^(Gamma_s)/(1-rho^2)
         d.Gamma_u <- ((Gamma_s * rho^(Gamma_s) * (1-rho^2)/rho) +
                      2 * rho * rho^(Gamma_s))/((1-rho^2)^2)
      } else {
         Gamma_u <- diag(TI)
         d.Gamma_u <- matrix(0, nrow=TI, ncol=TI)
         for (i in 2:TI) {
            d.Gamma_u[i-1, i] < -1
            d.Gamma_u[i, i-1] < -1 
         }
      }
      Gamma_v <- matrix(1,nrow=TI,ncol=TI)
      d.Gamma_v <- matrix(0,nrow=TI,ncol=TI)
   }
   u_corr <- diag(1, NV)
   if (NV > 1) {
      for (i in 2:NV) {
         for (j in 1:(i-1)){ 
            u_corr[i, j] <- u_corr[j, i] <- rho_u[((i-2)*(i-1))/2 + j]
         }
      }
   }
   sqrt_u <- sqrt(delta[1:NV])
   u_base <- ((sqrt_u %*% t(sqrt_u)) * u_corr) %x% Gamma_u
   sqrt_v <- sqrt(delta[(NV+1):(2*NV)])
   v_base <- ((sqrt_v %*% t(sqrt_v)) * u_corr) %x% Gamma_v 
   V <- vcov_e + diag(M) %x% (u_base + v_base)
   V.inv <- try(solve(V, tol=tol))
   if(inherits(V.inv, "try-error")) 
      stop("Initial V matrix (incl. random effects) is singular")
   V.inv.X <- try(solve(V, X, tol=tol))
   if(inherits(V.inv.X, "try-error"))
      stop("Initial V-inverse X matrix is singular")
   Xt.V.inv.X <- t(X) %*% V.inv.X
   B_Est   <- try(solve(Xt.V.inv.X, t(X) %*% solve(V, y, tol=tol),  # (6.2.5)
                  tol=tol))
   if(inherits(B_Est, "try-error"))
      stop("Initial calculation of beta failed")
   res <- y - X %*% B_Est
   A <- orth_c(X)
   llike.const <- N * (.5*log(2*pi) - mean(log(y.rescale)))
   if (method == "REML") {
      y.star <- t(A) %*% y
      C.y.star <- t(A) %*% V %*% A
      llike.c.const <- .5 * (N - dim(X)[2]) * log(2*pi) - 
                    N * mean(log(y.rescale)) + sum(ncolx*log(y.rescale))
      llike.c <- -.5 *(determinant(C.y.star)$modulus[1] +  
               t(y.star) %*% solve(C.y.star, y.star)) - llike.c.const
   } else {
      llike <- -.5 *(determinant(V)$modulus[1] +  t(res) %*%  
                   solve(V, res, tol=tol)) - llike.const
   }   
   conv.count <- 0
   bad.inf.mat <- FALSE
   # iterations over iter
   for (iter in 1:maxiter) {
      delta.hist[, iter] <- delta
      if (method == "REML") {
         llikelihood.hist[iter] <- llike.c
      } else {
         llikelihood.hist[iter] <- llike
      }        
      # compute the current information matrix corresponding to delta   
      V.inv.res <- try(solve(V, res, tol=tol))
      if(inherits(V.inv.res, "try-error")) {
         V.inv.res <- NA
         break
      }
      Xt.V.inv.X.inv <- try(solve(Xt.V.inv.X, tol=tol))
      if(inherits(Xt.V.inv.X.inv, "try-error")) {
         Xt.V.inv.X.inv <- NA
         break
      }
      P <- V.inv - V.inv.X %*% solve(Xt.V.inv.X, t(V.inv.X),   # (p. 101)
             tol=tol)
      s <- rep(0, times=len.delta)
      inf.mat <- matrix(0, nrow=len.delta, ncol=len.delta)

      V.j <- array(0 ,dim=c(TI*NV, TI*NV, len.delta))
      if (method == "REML") {   
         V.j.m <- array(0, dim=c(N, N, len.delta))
      } else {
         V.j.m <- array(0, dim=c(TI*NV, TI*NV, len.delta))
      }   
      for (i in 1:NV) { 
         c.m <- matrix(0, nrow=NV, ncol=NV)
         c.m[i, i] <- 1
         if(delta[i] > 0) {
            seq.j <- c(1:NV)[-i]
            for (j in seq.j) {
               if (j < i) {
                  c.m[i, j] <- rho_u[((i-2)*(i-1))/2+j]*
                              sqrt(delta[j]/delta[i])/2
               } else {
                  c.m[i, j] <- rho_u[((j-2)*(j-1))/2+i]*
                              sqrt(delta[j]/delta[i])/2
               }
               c.m[j, i] <- c.m[i, j]
            }
         }  
         V.j[, , i] <- c.m %x% Gamma_u
      }
     
      for (i in 1:NV) { 
         c.m <- matrix(0, nrow=NV, ncol=NV)
         c.m[i, i] <- 1
         if(delta[i+NV] > 0) {
            seq.j <- c(1:NV)[-i]
            for (j in seq.j) {
               if (j < i) {
                  c.m[i, j] <- rho_u[((i-2)*(i-1))/2+j]*
                              sqrt(delta[j+NV]/delta[i+NV])/2
               } else {
                  c.m[i, j] <- rho_u[((j-2)*(j-1))/2+i]*
                              sqrt(delta[j+NV]/delta[i+NV])/2
               }
               c.m[j, i] <- c.m[i, j]
            }
         } 
         V.j[, , i+NV] <- c.m %x% Gamma_v
      }

      V.j[, , index.rho]  <- 
                   ((sqrt_u %*% t(sqrt_u)) * u_corr) %x% d.Gamma_u +
                   ((sqrt_v %*% t(sqrt_v)) * u_corr) %x% d.Gamma_v
      if(NV > 1) {
         for (i in 2:NV) {
            for (j in 1:(i-1)) { 
               k <- ((i-2)*(i-1))/2+j
               c.m <- matrix(0, nrow=NV, ncol=NV)
               c.m[i, j] <- c.m[j, i] <- 1
               V.j[, , index.rho+k] <- 
                   ((sqrt_u %*% t(sqrt_u)) * c.m) %x% Gamma_u +
                   ((sqrt_v %*% t(sqrt_v)) * c.m) %x% Gamma_v
            }
         }
      }
      if (method == "REML") {     
         for (i in 1:len.delta) {
            for (m in 1:M) {
               r1 <- (m-1)*TI*NV + 1
               r2 <- m*TI*NV
               for (mp in 1:M) {
                  r1p <- (mp-1)*TI*NV + 1
                  r2p <- mp*TI*NV
                  V.j.m[r1:r2, r1p:r2p, i] <- P[r1:r2, r1p:r2p] %*% 
                                  as.matrix(V.j[, , i]) 
               }
               s[i] <- s[i] + .5 *
                    t(V.inv.res[r1:r2]) %*% as.matrix(V.j[, , i]) %*% 
                       V.inv.res[r1:r2]
            }
            s[i] <- s[i] - .5* sum(diag(V.j.m[, , i])) 
         }
         for (i in 1:len.delta) {
            for (j in (1:i)) {
               inf.mat[i, j] <- .5* sum(sapply(1:N,
                                   FUN=function(k) {
                                       t(as.vector(V.j.m[k, , i])) %*% 
                                         as.vector(V.j.m[, k, j])} ))  # (6.2.19)
               if(i != j) inf.mat[j, i] <- inf.mat[i, j]
            }
         }
      } else {
         for (m in 1:M) {
            r1 <- (m-1)*TI*NV+1
            r2 <- m*TI*NV
            V.inv <- solve(V[r1:r2, r1:r2], tol=tol)
            V.inv.res <- V.inv %*% res[r1:r2]
            for (i in 1:len.delta) {
               V.j.m[, , i] <- V.inv %*% as.matrix(V.j[, , i])
               s[i] <- s[i] -.5*sum(diag(as.matrix(V.j.m[, , i]))) + .5 *
                       t(V.inv.res) %*% as.matrix(V.j[, , i]) %*% V.inv.res 
            }
            for (i in 1:len.delta) {
               for (j in (1:i)) { 
                  inf.mat[i, j] <- inf.mat[i, j] +
                          .5* sum(diag(V.j.m[, , i] %*% V.j.m[, , j])) # (6.2.19)
               }
            }
         }
         for (i in 1:len.delta) for (j in (1:i)) { 
            if (i!=j) inf.mat[j, i] <- inf.mat[i, j]
         }
      }
      inf.mat.hist[, , iter] <- inf.mat
      s.hist[, iter] <- s
      if(llike.only) {
         res <- y - X %*% B_Est
         if (method == "REML") {       
            llike <- c(-.5 *(determinant(V)$modulus[1] +  t(res) %*%  
                       solve(V, res, tol=tol))) - llike.const
            if(NV > 1) {
               parm <- c(rho, delta[1:(2*NV)], rho_u, loglikelihood=llike,
                         constrained.ll=llike.c, num.iter=0)
            } else {
               parm <- c(rho, delta[1:(2*NV)], loglikelihood=llike,
                      constrained.ll=llike.c, num.iter=0)
            }
         } else {
            if(NV > 1) {
               parm <- c(rho, delta[1:(2*NV)], rho_u, loglikelihood=llike,
                         num.iter=0)
            } else {
               parm <- c(rho, delta[1:(2*NV)], loglikelihood=llike,
                         num.iter=0)
            }
         } 
         B_Est <- c(B_Est)
         names(B_Est) <- attr(X, "dimnames")[[2]]
         Xt.V.inv.X.inv <- try(solve(Xt.V.inv.X, tol=tol))
         if(inherits(Xt.V.inv.X.inv, "try-error")) {
            warning("Unable to invert X-transpose-V-inverse-X matrix", 
               immediate.=TRUE)
            Xt.V.inv.X.inv <- NA
         }
         if(is.null(rho.fixed)) {
	    test.inv <- try(solve(inf.mat, tol=tol))
	 } else {
	    test.inv <- try(solve(inf.mat[ix.start, ix.start], tol=tol))
	 }
         if(inherits(test.inv, "try-error")) 
            warning("Information matrix has become singular", immediate.=TRUE)
         if(use.y.rescale) {
            parm[2:(1+2*NV)] <- parm[2:(1+2*NV)]/(delta.adjust[1:(2*NV)]^2)
            B_Est <- B_Est/b.adjust
            delta[1:(2*NV)] <- delta[1:(2*NV)]/(delta.adjust[1:(2*NV)]^2)
            inf.mat <- delta.adjust^2 * t(inf.mat * delta.adjust^2)
            Xt.V.inv.X.inv <- (1/b.adjust) * t(Xt.V.inv.X.inv/b.adjust)
         }
         return(list(fit=fit, parm=parm, coef=B_Est, 
                     delta=delta, inf.mat=inf.mat, var.coef=Xt.V.inv.X.inv))
      }
      if(starting.delta & maxiter == 1) break
      if(iter > 1 & iter == maxiter) break
      if(iter > 5 | (iter > 1 & starting.delta) ) {
         if(max(abs(delta[1:(2*NV)] - last.delta[1:(2*NV)])) < 
            iter.tol*min.vcov &
            max(abs(delta[index.rho:len.delta] - 
                last.delta[index.rho:len.delta])) < iter.tol) {
            if(conv.count >= 1) {
               break 
            } else {
               conv.count <- conv.count + 1
            }
         } else {
            conv.count <- 0
         }
      }
      if(is.null(rho.fixed)) {
         inf.mat.inv <- try(solve(inf.mat, tol=tol))
      } else {
         inf.mat.inv <- try(solve(inf.mat[ix.start, ix.start], tol=tol))
      }
      if(inherits(inf.mat.inv, "try-error")) {
         bad.inf.mat <- TRUE
         print("Information matrix has become ill-conditioned")
         print(inf.mat)
         print("Current values of parameters (delta)")
         print(delta)
         print(paste("At iteration", iter))
         if(use.y.rescale)
            print("Note: Values reflect adjustment by y.rescale")
         break
      }
      last.delta <- delta
      if (method == "REML") {
         last.llike <- llike.c
      } else {
         last.llike <- llike
      }
      if(is.null(rho.fixed)) {
         delta.delta.base <- solve(inf.mat, s, tol=tol)
      } else {
         delta.delta.base <- rep(0, len.delta)
         delta.delta.base[ix.start] <- 
                solve(inf.mat[ix.start, ix.start], s[ix.start], tol=tol)
      }  
      if(iter < 5 & baby.steps) {
         delta.delta.base <- 2^(iter-5) * delta.delta.base 
      } else if(iter == 5 & baby.steps) {
         delta.delta.base <- .75*delta.delta.base 
      } 
      low.change <- TRUE
     
# iterations over adj.iter

      for(adj.iter in 1:20) {
         adj.ratio <- rep(1, times=len.delta)
         delta.delta <- delta.delta.base
         ix <- ix.start
         for(ix.index in 1:ix.len) { # iterations over ix.index
            if(NV > 1) for (i in c((index.rho+1):len.delta)) {
               if(i %in% ix) {
                  if (delta.delta[i] > 0) { 
                     adj.ratio[i] <- min(1,
                              (max.rho_u-last.delta[i])/delta.delta[i],
                               .05/delta.delta[i])
                  } else if (delta.delta[i] < 0 ) {
                    adj.ratio[i] <- min(1, 
                              (delta.min[i] - last.delta[i])/delta.delta[i],
                               -.05/delta.delta[i])
                  }                                   
               }
            }
            if (index.rho %in% ix) {
               if(!is.null(max.rho)) {
                  if (delta.delta[index.rho] > 0) {
                      adj.ratio[index.rho] <- min(1, 
                              (max.rho - last.delta[index.rho])/
                                        delta.delta[index.rho],
                              .05/delta.delta[index.rho])
                  } else if (delta.delta[index.rho] < 0 ) {
                      adj.ratio[index.rho] <- min(1, 
                              (delta.min[index.rho] - last.delta[index.rho])/
                                         delta.delta[index.rho],
                              -.05/delta.delta[index.rho])
                  }
               } else {
                  if (delta.delta[index.rho] > 0) {
                      adj.ratio[index.rho] <- min(1, 
                                    .05/delta.delta[index.rho])
                  } else if (delta.delta[index.rho] < 0 ) {
                      adj.ratio[index.rho] <- min(1, 
                          (delta.min[index.rho] - last.delta[index.rho])/
                                         delta.delta[index.rho],
                                    -.05/delta.delta[index.rho])
                  }
               }            
            }    
            for (i in 1:(2*NV)) {
               if (i %in% ix) {
                  if (delta.delta[i] < 0) {
                      adj.ratio[i] <- min(1,
                               (delta.min[i]-last.delta[i])/delta.delta[i])
                  }                    
               }
            }
            i.out <- which.min(adj.ratio[ix])

            if(adj.ratio[ix[i.out]] <= tol) {
               ix <- ix[-i.out]
               delta.delta[] <- 0
               inf.mat.inv <- try(solve(inf.mat[ix, ix]))
               if(inherits(inf.mat.inv, "try-error"))  next
               if(ix.index < ix.len) delta.delta[ix] <- (.5^(adj.iter-1)) *
                                solve(inf.mat[ix, ix], s[ix], tol=tol)
            } else if (length(ix) == 1 | adj.ratio[ix[i.out]] > .5^(iter/4) | 
                      iter%%2 == 0) {
               if(any(delta.delta > iter.tol)) low.change <- FALSE
               delta.delta <- adj.ratio[ix[i.out]] * delta.delta
               adj.factor.hist[iter] <- adj.ratio[ix[i.out]]
               break  
            } else {
               ix <- ix[-i.out]
               delta.delta[] <- 0
               inf.mat.inv <- try(solve(inf.mat[ix, ix]))
               if(inherits(inf.mat.inv, "try-error"))  next
               if(ix.index < ix.len) delta.delta[ix] <- (.5^(adj.iter-1))*
                                solve(inf.mat[ix, ix], s[ix], tol=tol)
            }
            
         } # end of loop over ix.index  
         if(length(ix)==0) {
            adj.factor.hist[iter] <- 1
            ix.hist[, iter] <- 0
            adj.hist[iter] <- adj.iter
            delta.delta.base <- .5*delta.delta.base
            next # next adj.iter
         }
         ix.hist[, iter] <- 0
         ix.hist[ix, iter] <- ix
         adj.hist[iter] <- adj.iter
         delta <- last.delta + dampening * delta.delta

         rho <- delta[index.rho]
         rho_u <- delta[(index.rho+1):len.delta]
         if (model == "dyn") {
            Gamma.list <- Gamma_u_f(rho, TI)
            Gamma_u <- Gamma.list[[1]]
            d.Gamma_u <- Gamma.list[[2]]
            Gamma.list <- Gamma_v_f(rho, TI)
            Gamma_v <- Gamma.list[[1]]
            d.Gamma_v <- Gamma.list[[2]]
         } else {
            Gamma_b <- matrix(rep(1:TI,times=TI),nrow=TI, ncol=TI)
            Gamma_s <- abs(Gamma_b-t(Gamma_b))
            if(rho > 0) {
               Gamma_u <- rho^(Gamma_s)/(1-rho^2)
               d.Gamma_u <- ((Gamma_s * rho^(Gamma_s) * (1-rho^2)/rho) +
                      2*rho* rho^(Gamma_s))/((1-rho^2)^2)
            } else {
               Gamma_u <- diag(TI)
               d.Gamma_u <- matrix(0,nrow=TI,ncol=TI)
               for (i in 2:TI) {
                  d.Gamma_u[i-1,i]<-1
                  d.Gamma_u[i,i-1]<-1
               }
            }
         }
         u_corr <- diag(1,NV)
         if (NV > 1) for (i in 2:NV) for (j in 1:(i-1)){ 
            u_corr[i,j] <- u_corr[j,i] <- rho_u[((i-2)*(i-1))/2+j]
            }
         sqrt_u <- sqrt(delta[1:NV])
         u_base <- ((sqrt_u %*% t(sqrt_u)) * u_corr) %x% Gamma_u
         sqrt_v <- sqrt(delta[(NV+1):(2*NV)])
         v_base <- ((sqrt_v %*% t(sqrt_v)) * u_corr) %x% Gamma_v
        
         V <- vcov_e + diag(M) %x% (u_base + v_base)
         V.inv <- try(solve(V, tol=tol))
         if(inherits(V.inv, "try-error")) {
            V.inv <- NA
            next
         }
         V.inv.X <- solve(V, X, tol=tol)
         Xt.V.inv.X <- t(X) %*% V.inv.X         
         B_Est <- try(solve(Xt.V.inv.X, t(X) %*% solve(V, y, tol=tol), # (6.2.5)
                         tol=tol))
         if(inherits(B_Est, "try-error")) {
            B_Est <- NA
            next
         }          
         res <- y - X %*% B_Est
         if (method == "REML") {
            C.y.star <- t(A) %*% V %*% A
            llike.c <- -.5 *(determinant(C.y.star)$modulus[1] +  
                    t(y.star) %*% solve(C.y.star,y.star)) - llike.c.const
            if (llike.c > last.llike) break
         } else {
            llike <- -.5 *(determinant(V)$modulus[1] +  t(res) %*%  
                    solve(V,res,tol=tol)) -  llike.const
            if (llike > last.llike) break
         }
         delta.delta.base <- .5 * delta.delta.base
      }
      if(any(is.na(V.inv))) {
         if(adj.iter == 20) {
            warning("V has become singular", immediate.=TRUE)
         } else {
            warning("During iteration, V became singular")
         }
         warning.hist[3, iter] <- 1
         break
      }
      if(any(is.na(B_Est))) {
         if(adj.iter == 20) {
            warning("Beta can no longer be estimated", immediate.=TRUE)
         } else {
            warning("During iteration, B could not be estimated")
         }
         warning.hist[4, iter] <- 1
         break
      }
      if(adj.iter == 20) {
         if (method == "REML") {
            if(last.llike - tol > llike.c & !low.change) {
               warning.hist[1, iter] <- 1
               warning(paste("At iteration", iter,
                   "no detected gain in restricted likelihood"))
            }
         } else {
            if(last.llike - tol > llike & !low.change) {
               warning.hist[1, iter] <- 1
               warning(paste("At iteration", iter,
                   "no detected gain in likelihood"))
            }
         }
      }
   }
   if(any(is.na(V.inv)) | any(is.na(B_Est)) | bad.inf.mat){
      length(llikelihood.hist) <- iter
      delta.hist <- delta.hist[, 1:iter, drop=FALSE]
      length(adj.hist) <- iter
      length(adj.factor.hist) <- iter
      inf.mat.hist <- inf.mat.hist[, , 1:iter, drop=FALSE]
      s.hist <- s.hist[, 1:iter, drop=FALSE]
      ix.hist <- ix.hist[, 1:iter, drop=FALSE]
      adj.factor.hist < adj.factor.hist[1:iter]
      warning.hist <- warning.hist[, 1:iter, drop=FALSE]
      if(use.y.rescale) {
         delta[1:(2*NV)] <- delta[1:(2*NV)]/(delta.adjust[1:(2*NV)]^2)
         inf.mat <- delta.adjust^2 * t(inf.mat * delta.adjust^2)
         delta.hist[1:(2*NV), ] <- delta.hist[1:(2*NV), ]/
                                  (delta.adjust[1:(2*NV)]^2)
         for (i in 1:iter) {
            inf.mat.hist[, , i] <- delta.adjust^2 * 
                                t(inf.mat.hist[, , i] * delta.adjust^2)
         }   
#         s.hist[1:(2*NV), ] <- delta.adjust[1:(2*NV)] * s.hist[1:(2*NV), ] 
      }
      return(list(fit=fit, delta=delta, inf.mat=inf.mat,
             delta.hist=delta.hist,
             llikelihood.hist=llikelihood.hist,
             adj.hist=adj.hist,
             inf.mat.hist=inf.mat.hist,
             s.hist=s.hist,
             ix.hist=ix.hist,
             adj.factor.hist=adj.factor.hist,
             warning.hist=warning.hist))
   }

   if(iter == maxiter & maxiter > 1) {
      warning.hist[2, iter] <- 2
      warning("maxiter=", maxiter, " reached, check convergence",
         immediate. = TRUE)
   } else if (conv.count >= 1) {
      fit$convergence <- TRUE
   }
   if (method == "REML") 
      llike <- c(-.5 *(determinant(V)$modulus[1] +  t(res) %*%  
                       solve(V, res, tol=tol))) - llike.const
   V_bar <- solve(inf.mat, tol=tol)
   if(use.y.rescale) {
      B_Est <- B_Est/b.adjust
      delta[1:(2*NV)] <- delta[1:(2*NV)]/(delta.adjust[1:(2*NV)]^2)
      inf.mat <- delta.adjust^2 * t(inf.mat * delta.adjust^2)
      V_var <- (1/delta.adjust^2) * t(V_bar * (1/delta.adjust^2))
      Xt.V.inv.X.inv <- (1/b.adjust) * t(Xt.V.inv.X.inv/b.adjust)
      vcov_e <- (1/y.adjust) * t(vcov_e/y.adjust)
      y <- y/y.adjust
      V <- (1/y.adjust) * t(V/y.adjust)
   }

   sqrt_u <- sqrt(delta[1:NV])
   u_base <- ((sqrt_u %*% t(sqrt_u)) * u_corr) %x% Gamma_u
   sqrt_v <- sqrt(delta[(NV+1):(2*NV)])
   v_base <- ((sqrt_v %*% t(sqrt_v)) * u_corr) %x% Gamma_v
   M_term_base <- u_base + v_base
   for (m in 1:M){
      r1 <- NV*(m-1)*TI + 1
      r2 <- NV*m*TI
      Xi <- X[r1:r2, ]
      for (comp in 1:n.contr) {
         mi <- mi.mat[, comp]
         li <- mi %*% Xi
         M_term.t <- t(mi) %*% M_term_base
         contrast.fixed.est[m, comp] <- li %*% B_Est 
         contrast.est[m, comp] <- contrast.fixed.est[m, comp] +
             M_term.t %*% 
             solve(V[r1:r2, r1:r2], (y[r1:r2] - X[r1:r2,]%*%B_Est), tol=tol)
         contrast.wt1[m, comp] <- t(mi) %*% 
             solve(V[r1:r2, r1:r2], t(M_term.t), tol=tol) / sum(mi*mi)
         contrast.wt2[m, comp] <- contrast.wt1[m,comp] + t(mi) %*% 
             (diag(NV*TI) - solve(V[r1:r2,r1:r2], M_term_base, tol=tol)) %*%
             (Xi %*% Xt.V.inv.X.inv %*% t(Xi) %*% 
             solve(V[r1:r2,r1:r2], mi, tol=tol))/sum(mi*mi)
      }
   }
 

# begin MSE calculation
   Gi <- u_base + v_base
   d.bi <- matrix(0, nrow=NV*TI, ncol=len.delta)
   d.V  <- array(0, dim=c(NV*TI, NV*TI, len.delta))
   g1 <- g2 <- g3 <- fixed.var <- matrix(0, nrow=M, ncol=n.contr)

   for (i in 1:NV) { 
      c.m <- matrix(0, nrow=NV, ncol=NV)
      c.m[i, i] <- 1
      if(delta[i] > 0) {
         seq.j <- c(1:NV)[-i]
         for (j in seq.j) {
            if (j < i) {
               c.m[i, j] <- rho_u[((i-2)*(i-1))/2+j]*
                           sqrt(delta[j]/delta[i])/2
            } else {
               c.m[i, j] <- rho_u[((j-2)*(j-1))/2+i]*
                           sqrt(delta[j]/delta[i])/2
            }
            c.m[j, i] <- c.m[i, j]
         }
      }  
      d.V[, , i]<- c.m %x% Gamma_u
   }

   for (i in 1:NV) { 
      c.m <- matrix(0, nrow=NV, ncol=NV)
      c.m[i, i] <- 1
      if(delta[i] > 0) {
         seq.j <- c(1:NV)[-i]
         for (j in seq.j) {
            if (j < i) {
               c.m[i, j] <- rho_u[((i-2)*(i-1))/2+j]*
                           sqrt(delta[j+NV]/delta[i+NV])/2
            } else {
               c.m[i, j] <- rho_u[((j-2)*(j-1))/2+i]*
                           sqrt(delta[j+NV]/delta[i+NV])/2
            }
            c.m[j, i] <- c.m[i, j]
         }
      }  
      d.V[, , i+NV] <- c.m %x% Gamma_v
   }

   d.V[, , index.rho] <- ((sqrt_u%*%t(sqrt_u)) * u_corr) %x% d.Gamma_u +
                         ((sqrt_v%*%t(sqrt_v)) * u_corr) %x% d.Gamma_v
             
   if(NV > 1) for (i in 2:NV) for (j in 1:(i-1)){ 
      k <- ((i-2)*(i-1))/2 + j
      c.m <- matrix(0, nrow=NV, ncol=NV)
      c.m[i, j] <- c.m[j, i] <- 1
      d.V[, , index.rho+k] <- 
                 ((sqrt_u%*%t(sqrt_u)) * c.m) %x% Gamma_u +
                 ((sqrt_v%*%t(sqrt_v)) * c.m) %x% Gamma_v
      }
 
   for (m in 1:M){
      r1 <- NV*(m-1)*TI + 1
      r2 <- NV*m*TI
      Xi <- X[r1:r2, ]
      Vi <- V[r1:r2, r1:r2]
      vcov_ei <- vcov_e[r1:r2, r1:r2]
      Vi.inv  <- solve(Vi, tol=tol)
      for (comp in c(1:n.contr)) {
         mi <- mi.mat[, comp]
         li <- mi %*% Xi
         bi <- t(mi) %*% Gi %*% Vi.inv                     # (6.2.10)
         di <- li - bi %*% Xi                              # (6.2.10)
         for (i in 1:len.delta) {
            d.bi[, i] <- t(t(mi) %*% as.matrix(d.V[, , i]) %*% Vi.inv - 
                         t(mi)  %*% Gi %*% Vi.inv %*% 
                         as.matrix(d.V[, , i]) %*% Vi.inv)
         }
         g1[m, comp] <- t(mi) %*% (Gi - Gi %*% Vi.inv %*% Gi) %*% mi
         g2[m, comp] <- di %*% Xt.V.inv.X.inv %*% t(di)
         g3[m, comp] <- sum(diag((t(d.bi) %*% Vi %*% d.bi) %*% V_bar))
         fixed.var[m, comp] <- li %*% Xt.V.inv.X.inv %*% t(li)
      }
   }
  
      
   if(!is.null(y.include)) {
      y.out <- which(y.include != 1) 
      MM <- length(y.out)
      if(MM > 0){
         merge.c <- function(y1, y2) {
            y <- rep(0, M.save)         
            y[y.in] <- y1
            y[y.out] <- y2
            return(y)
         }
         merge.m <- function(y1, y2) {
            y <- matrix(0, nrow=M.save, ncol=dim(y1)[2])
            for (i in 1:dim(y1)[2]) {
               y[, i] <- merge.c(as.vector(y1[, i]),
                                 as.vector(y2[, i]))
            }
            return(y)
         }
         index.non.obs <- rep(0, NV*MM*TI)
         g1.non <- g2.non <- zero.non <- matrix(0, nrow=MM, ncol=n.contr)
         contrast.est.non <- contrast.fixed.non <- 
                             matrix(0, nrow=MM, ncol=n.contr)
         for (m in 1:MM) {
            index.non.obs[(NV*TI*(m-1)+1):(NV*TI*m)] <-
                        c((NV*TI*(y.out[m]-1)+1):(NV*TI*y.out[m]))
            Xi <- as.matrix(X.save[(NV*TI*(y.out[m]-1) + 1):(NV*TI*y.out[m]),])
            for (comp in c(1:n.contr)) {
               mi <- mi.mat[, comp]
               li <- mi %*% Xi
               g1.non[m, comp] <- t(mi) %*% Gi %*%mi
               g2.non[m, comp] <- li %*% Xt.V.inv.X.inv %*% t(li) 
               contrast.est.non[m, comp] <- li %*% B_Est
            } 
         }
         contrast.est <- merge.m(contrast.est, contrast.est.non)
         contrast.fixed.est <- merge.m(contrast.fixed.est, contrast.est.non)
         contrast.wt1 <- merge.m(contrast.wt1, zero.non)
         contrast.wt2 <- merge.m(contrast.wt2, zero.non)
         g1 <- merge.m(g1, g1.non)
         g2 <- merge.m(g2, g2.non)
         g3 <- merge.m(g3, zero.non)
         fixed.var <- merge.m(fixed.var, g2.non)
         M <- M.save
      }
   } 
  
   if(n.contr > NV*TI) {
      contrast.g1  <- as.matrix(g1[, (NV*TI+1):n.contr])
      contrast.g2  <- as.matrix(g2[, (NV*TI+1):n.contr])
      contrast.g3  <- as.matrix(g3[, (NV*TI+1):n.contr])
      contrast.mse <- contrast.g1 + contrast.g2 + 2*contrast.g3
      contrast.fixed.var <- as.matrix(fixed.var[, (NV*TI+1):n.contr])
   }
   if(NV > 1) {
      eblup <- eblup.g1 <- eblup.g2 <- eblup.g3 <- eblup.wt1 <-
               eblup.wt2 <- est.fixed <- est.fixed.var <-
               matrix(0, nrow=M*TI, ncol=NV)
      for (nv in 1:NV) {
         for (t in 1:TI) {
            eblup[TI*(0:(M-1))+t, nv] <-
                      contrast.est[1:M, TI*(nv-1) + t]
            est.fixed[TI*(0:(M-1))+t, nv] <-
                      contrast.fixed.est[1:M, TI*(nv-1) + t]
            eblup.wt1[TI*(0:(M-1))+t, nv] <-
                      contrast.wt1[1:M, TI*(nv-1) + t]
            eblup.wt2[TI*(0:(M-1))+t, nv] <-
                      contrast.wt2[1:M, TI*(nv-1) + t]
            eblup.g1[TI*(0:(M-1))+t, nv] <-
                      g1[1:M, TI*(nv-1) + t]
            eblup.g2[TI*(0:(M-1))+t, nv] <-
                      g2[1:M, TI*(nv-1) + t]
            eblup.g3[TI*(0:(M-1))+t, nv] <-
                      g3[1:M,TI*(nv-1) + t]
            est.fixed.var[TI*(0:(M-1))+t, nv] <-
                      fixed.var[1:M, TI*(nv-1) + t]
         }
      }
   } else {
      eblup <- eblup.g1 <- eblup.g2 <- eblup.g3 <- eblup.wt1 <-
               eblup.wt2 <- est.fixed <- est.fixed.var <- rep(0, times=M*TI)
      for (t in 1:TI) {
         eblup[TI*(0:(M-1))+t] <- contrast.est[1:M, t]
         est.fixed[TI*(0:(M-1))+t] <- contrast.fixed.est[1:M, t]
         eblup.wt1[TI*(0:(M-1))+t] <- contrast.wt1[1:M, t]
         eblup.wt2[TI*(0:(M-1))+t] <- contrast.wt2[1:M, t]
         eblup.g1[TI*(0:(M-1))+t] <- g1[1:M, t]
         eblup.g2[TI*(0:(M-1))+t] <- g2[1:M, t]
         eblup.g3[TI*(0:(M-1))+t] <- g3[1:M, t]
         est.fixed.var[TI*(0:(M-1))+t] <- fixed.var[1:M, t]
      }
   }
   eblup.mse <- eblup.g1 + eblup.g2 + 2*eblup.g3
  
   if(n.contr > NV*TI) {
      contrast.fixed.est <- as.matrix(contrast.fixed.est[, (NV*TI+1):n.contr])
      contrast.est <- as.matrix(contrast.est[, (NV*TI+1):n.contr])
      contrast.wt1 <- as.matrix(contrast.wt1[, (NV*TI+1):n.contr]) 
      contrast.wt2 <- as.matrix(contrast.wt2[, (NV*TI+1):n.contr])
   } else {
      contrast.fixed.est <- contrast.est <- contrast.wt1 <- contrast.wt2 <- 
          NULL
   }
   if (method == "REML") {
      if(NV > 1) {
         parm <- c(rho, delta[1:(2*NV)], rho_u, loglikelihood=llike,
                   constrained.ll=llike.c, num.iter=iter)
      } else {
         parm <- c(rho, delta[1:(2*NV)], loglikelihood=llike,
                   constrained.ll=llike.c, num.iter=iter)
      }
   } else {
      if(NV > 1) {
         parm <- c(rho, delta[1:(2*NV)], rho_u, loglikelihood=llike,
                   num.iter=iter)
      } else {
         parm <- c(rho, delta[1:(2*NV)], loglikelihood=llike,
                   num.iter=iter)
      }
   }
   B_Est <- c(B_Est)
   names(B_Est) <- attr(X, "dimnames")[[2]]
   keep.history <- FALSE
   if (is.null(iter.history)) {
      if(!fit$convergence & maxiter > 1) 
         keep.history <- TRUE
   } else {
      if(iter.history)
         keep.history <- TRUE
   }      
   if (keep.history) {
      length(llikelihood.hist) <- iter
      delta.hist <- delta.hist[, 1:iter, drop=FALSE]
      length(adj.hist) <- iter
      inf.mat.hist <- inf.mat.hist[, , 1:iter, drop=FALSE]
      s.hist <- s.hist[, 1:iter, drop=FALSE]
      ix.hist <- ix.hist[, 1:iter, drop=FALSE]
      adj.factor.hist <- adj.factor.hist[1:iter]
      warning.hist <- warning.hist[, 1:iter, drop=FALSE]
      if(use.y.rescale) {
         delta.hist[1:(2*NV), ] <- delta.hist[1:(2*NV), ]/
                                  (delta.adjust[1:(2*NV)]^2)
         for (i in 1:iter) {
            inf.mat.hist[, , i] <- delta.adjust^2 * 
                                t(inf.mat.hist[, , i] * delta.adjust^2)
         }
      }
      return(list(eblup=eblup, fit=fit, parm=parm, coef=B_Est,  ids=ids,  
                delta=delta, eblup.mse=eblup.mse, eblup.g1=eblup.g1,   
                eblup.g2=eblup.g2, eblup.g3=eblup.g3, est.fixed=est.fixed, 
                est.fixed.var=est.fixed.var, eblup.wt1=eblup.wt1, 
                eblup.wt2=eblup.wt2, contrast.est=contrast.est, 
                contrast.mse=contrast.mse, contrast.g1=contrast.g1, 
                contrast.g2=contrast.g2, contrast.g3=contrast.g3,
                contrast.fixed.est=contrast.fixed.est,
                contrast.fixed.var=contrast.fixed.var,
                contrast.wt1=contrast.wt1, contrast.wt2=contrast.wt2,
                inf.mat=inf.mat, var.coef=Xt.V.inv.X.inv,
                delta.hist=delta.hist,
                llikelihood.hist=llikelihood.hist,
                adj.hist=adj.hist,
                inf.mat.hist=inf.mat.hist,
                s.hist=s.hist,
                ix.hist=ix.hist,
                adj.factor.hist=adj.factor.hist,
                warning.hist=warning.hist))
   } else {
    
      return(list(eblup=eblup, fit=fit, parm=parm, coef=B_Est, ids=ids, 
                delta=delta, eblup.mse=eblup.mse, eblup.g1=eblup.g1,   
                eblup.g2=eblup.g2, eblup.g3=eblup.g3, est.fixed=est.fixed, 
                est.fixed.var=est.fixed.var, eblup.wt1=eblup.wt1, 
                eblup.wt2=eblup.wt2, contrast.est=contrast.est, 
                contrast.mse=contrast.mse, contrast.g1=contrast.g1, 
                contrast.g2=contrast.g2, contrast.g3=contrast.g3,
                contrast.fixed.est=contrast.fixed.est,
                contrast.fixed.var=contrast.fixed.var,
                contrast.wt1=contrast.wt1, contrast.wt2=contrast.wt2,
                inf.mat=inf.mat, var.coef=Xt.V.inv.X.inv ))
   }
}
